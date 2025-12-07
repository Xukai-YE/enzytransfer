#!/usr/bin/env python3
"""
序列对比与合并工具 - 并行优化版本
使用多进程加速文件处理
"""

import yaml
import os
import sys
import argparse
from pathlib import Path
from collections import defaultdict
import difflib
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import hashlib

# 尝试使用更快的CLoader
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
    print("⚠️  提示: 安装libyaml可以加速YAML处理 (pip install pyyaml --global-option='--with-libyaml')")

def compute_sequence_hash(sequence):
    """计算序列的hash值，用于快速比较"""
    if not sequence:
        return None
    return hashlib.md5(sequence.encode()).hexdigest()

def extract_info_from_yaml(yaml_path):
    """从YAML文件中提取关键信息（优化版）"""
    try:
        with open(yaml_path, 'r', encoding='utf-8') as f:
            # 使用CLoader加速
            data = yaml.load(f, Loader=Loader)
        
        proteins = data.get('proteins', [])
        if not proteins:
            return None
        
        protein = proteins[0]
        
        wildtype_seq = protein.get('wildtype_sequence')
        mutant_seq = protein.get('mutant_sequence')
        
        info = {
            'file_path': str(yaml_path),
            'file_name': yaml_path.name,
            'uniprotid': protein.get('uniprotid'),
            'variant_type': protein.get('variant_type'),
            'variant_description': protein.get('variant_description'),
            'wildtype_sequence': wildtype_seq,
            'mutant_sequence': mutant_seq,
            'wildtype_hash': compute_sequence_hash(wildtype_seq),  # 用于快速比较
            'mutant_hash': compute_sequence_hash(mutant_seq),
            'organism': protein.get('organism'),
            'ecnumber': protein.get('ecnumber'),
            'full_data': data
        }
        
        return info
    
    except Exception as e:
        print(f"⚠️  读取文件出错 {yaml_path}: {e}")
        return None

def scan_yaml_files(directories):
    """扫描指定目录下的所有YAML文件"""
    yaml_files = []
    
    for directory in directories:
        dir_path = Path(directory)
        if not dir_path.exists():
            print(f"⚠️  目录不存在: {directory}")
            continue
        
        for yaml_file in dir_path.rglob('*.yaml'):
            yaml_files.append(yaml_file)
        for yaml_file in dir_path.rglob('*.yml'):
            yaml_files.append(yaml_file)
    
    return yaml_files

def process_files_parallel(yaml_files, max_workers=None):
    """
    并行处理YAML文件
    max_workers: 工作进程数，默认为CPU核心数
    """
    if max_workers is None:
        max_workers = min(cpu_count(), len(yaml_files))
    
    yaml_infos = []
    
    print(f"使用 {max_workers} 个进程并行处理...")
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_file = {
            executor.submit(extract_info_from_yaml, yaml_file): yaml_file 
            for yaml_file in yaml_files
        }
        
        # 收集结果（带进度显示）
        completed = 0
        total = len(yaml_files)
        
        for future in as_completed(future_to_file):
            completed += 1
            if completed % 10 == 0 or completed == total:
                print(f"  处理进度: {completed}/{total} ({completed*100//total}%)", end='\r')
            
            try:
                info = future.result()
                if info:
                    yaml_infos.append(info)
            except Exception as e:
                yaml_file = future_to_file[future]
                print(f"\n⚠️  处理文件失败 {yaml_file}: {e}")
        
        print()  # 换行
    
    return yaml_infos

def group_by_uniprot_mutation(yaml_infos):
    """按照uniprotid和mutation分组"""
    groups = defaultdict(list)
    
    for info in yaml_infos:
        if info and info['uniprotid']:
            mutation = info['variant_description'] if info['variant_description'] else 'wildtype'
            key = (info['uniprotid'], mutation)
            groups[key].append(info)
    
    return groups

def compare_sequences_fast(info_list):
    """
    快速序列比较（使用hash）
    返回: (是否一致, 消息)
    """
    if len(info_list) < 2:
        return True, "Only one file"
    
    # 先用hash快速比较
    wildtype_hashes = [f['wildtype_hash'] for f in info_list if f['wildtype_hash']]
    if len(wildtype_hashes) > 1:
        if len(set(wildtype_hashes)) > 1:
            return False, "Wildtype sequences differ"
    
    mutant_hashes = [f['mutant_hash'] for f in info_list if f['mutant_hash']]
    if len(mutant_hashes) > 1:
        if len(set(mutant_hashes)) > 1:
            return False, "Mutant sequences differ"
    
    return True, "All sequences are identical"

def merge_yaml_files(file_infos, output_dir):
    """合并多个YAML文件为一个"""
    
    if not file_infos:
        return None
    
    base_data = file_infos[0]['full_data'].copy()
    base_info = file_infos[0]
    
    merged_data = {
        'version': '2.0',
        'created': datetime.now().isoformat(),
        'modified': None,
        'creators': [],
        'vessels': None,
        'proteins': base_data.get('proteins', [])[:1],
        'complexes': None,
        'small_molecules': [],
        'reactions': [],
        'measurements': [],
        'equations': None,
        'parameters': [],
        'buffer': {'buffer': None, 'buffer_concentration': None},
        'references': []
    }
    
    # 合并creators
    creators_set = set()
    for info in file_infos:
        data = info['full_data']
        for creator in data.get('creators', []):
            creator_tuple = (
                creator.get('given_name', ''),
                creator.get('family_name', ''),
                creator.get('mail', '')
            )
            creators_set.add(creator_tuple)
    
    merged_data['creators'] = [
        {'given_name': c[0], 'family_name': c[1], 'mail': c[2]}
        for c in sorted(creators_set)
    ]
    
    # 合并small_molecules
    molecules_dict = {}
    for info in file_infos:
        data = info['full_data']
        for mol in data.get('small_molecules', []):
            mol_id = mol.get('id')
            if mol_id and mol_id not in molecules_dict:
                molecules_dict[mol_id] = mol
            elif mol_id and mol.get('name'):
                if not molecules_dict[mol_id].get('name'):
                    molecules_dict[mol_id] = mol
    
    merged_data['small_molecules'] = list(molecules_dict.values())
    
    # 合并reactions, measurements, parameters
    for info in file_infos:
        data = info['full_data']
        merged_data['reactions'].extend(data.get('reactions', []))
        merged_data['measurements'].extend(data.get('measurements', []))
        merged_data['parameters'].extend(data.get('parameters', []))
    
    # 合并references
    references_dict = {}
    for info in file_infos:
        data = info['full_data']
        for ref in data.get('references', []):
            key = ref.get('pubmed_id') or ref.get('title') or str(len(references_dict))
            if key not in references_dict:
                references_dict[key] = ref
    
    merged_data['references'] = list(references_dict.values())
    
    # 合并buffer
    buffer_info = []
    for info in file_infos:
        data = info['full_data']
        buffer_data = data.get('buffer', {})
        if buffer_data and buffer_data.get('buffer'):
            buffer_info.append(buffer_data.get('buffer'))
    
    if buffer_info:
        merged_data['buffer']['buffer'] = '; '.join(set(buffer_info))
    
    # 生成文件名
    uniprotid = base_info['uniprotid']
    mutation = base_info['variant_description'] or 'wildtype'
    safe_organism = (base_info['organism'] or 'Unknown').replace(' ', '_').replace('/', '_')
    safe_mutation = mutation.replace(' ', '_').replace('/', '_')
    
    merged_filename = f"MERGED_{uniprotid}_{safe_mutation}.yaml"
    
    merged_data['name'] = f"Merged entry for {uniprotid} {mutation}"
    merged_data['description'] = f"Merged from {len(file_infos)} files: " + ", ".join([f['file_name'] for f in file_infos])
    
    # 写入文件（使用CDumper加速）
    output_path = os.path.join(output_dir, merged_filename)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        yaml.dump(merged_data, f, Dumper=Dumper, default_flow_style=False, 
                  allow_unicode=True, sort_keys=False)
    
    return output_path

def main():
    parser = argparse.ArgumentParser(
        description='序列对比与合并工具 - 并行优化版',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python merge_sequences_fast.py /data/sabio /data/brenda -o merged
  python merge_sequences_fast.py /data/db1 /data/db2 -j 8  # 使用8个进程
        """
    )
    
    parser.add_argument('directories', nargs='+', 
                       help='要扫描的目录路径(至少需要2个目录)')
    parser.add_argument('-o', '--output-dir', default='merged_output', 
                       help='合并后的YAML文件输出目录 (默认: merged_output)')
    parser.add_argument('-j', '--jobs', type=int, default=None,
                       help='并行进程数 (默认: CPU核心数)')
    parser.add_argument('--allow-single-dir', action='store_true',
                       help='允许只扫描单个目录（用于测试）')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='显示详细输出')
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("序列对比与合并工具 - 并行优化版 🚀")
    print("=" * 80)
    print()
    
    # 检查目录数量
    if len(args.directories) < 2 and not args.allow_single_dir:
        print("❌ 错误: 合并工具需要至少2个目录作为输入!")
        print()
        print("合并工具的使用场景是整合来自不同数据库的数据，例如:")
        print("  python3 merge_sequences_fast.py /data/sabio_output /data/brenda_output")
        print()
        print("如果你只想检查单个目录，请使用对比工具:")
        print("  python3 compare_sequences_cli.py /data/your_folder")
        print()
        print("如果确实需要在单个目录中查找并合并重复项，请添加 --allow-single-dir 选项")
        sys.exit(1)
    
    if len(args.directories) == 1 and args.allow_single_dir:
        print("⚠️  警告: 只扫描单个目录（测试模式）")
        print()
    
    # 验证目录
    valid_dirs = []
    for directory in args.directories:
        if os.path.exists(directory):
            valid_dirs.append(directory)
            print(f"✓ 找到目录: {directory}")
        else:
            print(f"✗ 目录不存在: {directory}")
    
    if not valid_dirs:
        print("\n❌ 错误: 没有有效的目录!")
        sys.exit(1)
    
    # 创建输出目录
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n✓ 输出目录: {output_dir}")
    
    # 确定进程数
    max_workers = args.jobs or cpu_count()
    print(f"✓ 使用 {max_workers} 个CPU核心")
    
    print(f"\n{'='*80}")
    print("阶段 1/4: 扫描YAML文件")
    print('='*80)
    
    # 扫描文件
    yaml_files = scan_yaml_files(valid_dirs)
    print(f"找到 {len(yaml_files)} 个YAML文件")
    
    if not yaml_files:
        print("\n❌ 错误: 没有找到YAML文件!")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print("阶段 2/4: 并行提取文件信息")
    print('='*80)
    
    # 并行处理文件
    yaml_infos = process_files_parallel(yaml_files, max_workers=max_workers)
    
    print(f"✓ 成功提取 {len(yaml_infos)} 个文件的信息")
    
    if args.verbose:
        print("\n提取的文件列表:")
        for info in yaml_infos[:10]:  # 只显示前10个
            print(f"  {info['uniprotid']} | {info['variant_description']} | {info['file_name']}")
        if len(yaml_infos) > 10:
            print(f"  ... 还有 {len(yaml_infos) - 10} 个文件")
    
    print(f"\n{'='*80}")
    print("阶段 3/4: 分组和序列检查")
    print('='*80)
    
    # 按uniprotid和mutation分组
    groups = group_by_uniprot_mutation(yaml_infos)
    print(f"找到 {len(groups)} 个不同的 uniprotid+mutation 组合")
    
    # 找出有多个文件的组
    duplicate_groups = {k: v for k, v in groups.items() if len(v) > 1}
    
    print(f"发现 {len(duplicate_groups)} 个有重复的组合")
    
    if not duplicate_groups:
        print("\n没有发现重复的uniprotid+mutation组合,无需合并。")
        sys.exit(0)
    
    print(f"\n{'='*80}")
    print("阶段 4/4: 执行合并")
    print('='*80)
    
    # 对每个重复组进行处理
    merged_count = 0
    skipped_count = 0
    
    for idx, ((uniprotid, mutation), files) in enumerate(duplicate_groups.items(), 1):
        print(f"\n[{idx}/{len(duplicate_groups)}] {uniprotid} | {mutation} ({len(files)}个文件)")
        
        if args.verbose:
            for f in files:
                print(f"  - {f['file_name']}")
        
        # 快速序列检查（使用hash）
        is_identical, msg = compare_sequences_fast(files)
        
        if is_identical:
            print(f"  ✓ 序列一致，正在合并...")
            
            try:
                merged_file = merge_yaml_files(files, output_dir)
                if merged_file:
                    print(f"  ✓ 合并成功: {os.path.basename(merged_file)}")
                    merged_count += 1
                else:
                    print("  ✗ 合并失败")
                    skipped_count += 1
            except Exception as e:
                print(f"  ✗ 合并失败: {e}")
                skipped_count += 1
        else:
            print(f"  ⚠️  {msg}，跳过合并")
            skipped_count += 1
    
    # 总结
    print(f"\n{'=' * 80}")
    print("✅ 处理完成!")
    print("=" * 80)
    print(f"总共扫描: {len(yaml_files)} 个YAML文件")
    print(f"成功解析: {len(yaml_infos)} 个文件")
    print(f"发现重复组合: {len(duplicate_groups)} 个")
    print(f"成功合并: {merged_count} 个")
    print(f"跳过: {skipped_count} 个")
    print(f"输出目录: {output_dir}")
    print("=" * 80)

if __name__ == "__main__":
    main()
