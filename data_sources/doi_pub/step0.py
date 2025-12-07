"""
DOI到PMID转换工具
使用方法：在本地环境运行此脚本
需要安装：pip install pandas openpyxl requests
"""

import pandas as pd
import requests
import time
from typing import Optional

def get_pmid_from_doi(doi: str, email: str = "your_email@example.com") -> Optional[str]:
    """
    根据DOI获取PubMed ID (PMID)
    
    参数:
        doi: 文章的DOI (可以是URL格式或纯DOI)
        email: 你的邮箱地址(NCBI要求提供)
    
    返回:
        PMID字符串,如果未找到则返回None
    """
    # 清理DOI格式
    doi = str(doi).strip()
    if doi.lower() in ['nan', 'none', '']:
        return None
    
    # 从URL中提取纯DOI
    if 'doi.org/' in doi:
        doi = doi.split('doi.org/')[-1]
    
    # 使用NCBI E-utilities API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    params = {
        'db': 'pubmed',
        'term': f'{doi}[DOI]',
        'retmode': 'json',
        'email': email
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        
        # 检查是否找到结果
        if 'esearchresult' in data and 'idlist' in data['esearchresult']:
            id_list = data['esearchresult']['idlist']
            if id_list:
                return id_list[0]  # 返回第一个PMID
        
        return None
        
    except Exception as e:
        print(f"  ✗ 查询出错: {e}")
        return None


def process_excel_with_doi(input_file: str, output_file: str, 
                           doi_column: str = 'html_doi',
                           email: str = "your_email@example.com"):
    """
    处理包含DOI的Excel文件,添加PMID列
    
    参数:
        input_file: 输入Excel文件路径
        output_file: 输出Excel文件路径
        doi_column: DOI所在的列名
        email: 你的邮箱
    """
    print(f"正在读取文件: {input_file}")
    df = pd.read_excel(input_file)
    
    # 检查DOI列是否存在
    if doi_column not in df.columns:
        print(f"错误: 找不到列 '{doi_column}'")
        print(f"可用的列名: {df.columns.tolist()}")
        return
    
    # 获取所有唯一的DOI
    unique_dois = df[doi_column].dropna().unique()
    print(f"\n找到 {len(unique_dois)} 个不同的DOI")
    print(f"总共 {len(df)} 行数据\n")
    
    # 创建DOI到PMID的映射字典
    doi_to_pmid = {}
    
    # 逐个查询DOI
    for i, doi_url in enumerate(unique_dois, 1):
        # 提取纯DOI用于显示
        display_doi = doi_url.split('doi.org/')[-1] if 'doi.org/' in doi_url else doi_url
        print(f"[{i}/{len(unique_dois)}] 查询DOI: {display_doi}")
        
        pmid = get_pmid_from_doi(doi_url, email)
        doi_to_pmid[doi_url] = pmid
        
        if pmid:
            print(f"  ✓ 找到PMID: {pmid}")
        else:
            print(f"  ✗ 未找到PMID")
        
        # 避免超过API请求限制(每秒最多3个请求)
        if i < len(unique_dois):
            time.sleep(0.34)
    
    # 将PMID映射回原始DataFrame
    print("\n正在回填PMID到表格...")
    df['PMID'] = df[doi_column].map(doi_to_pmid)
    
    # 统计结果
    total_rows = len(df)
    found_pmid = df['PMID'].notna().sum()
    missing_pmid = df['PMID'].isna().sum()
    
    print(f"\n处理完成!")
    print(f"总行数: {total_rows}")
    print(f"找到PMID: {found_pmid} ({found_pmid/total_rows*100:.1f}%)")
    print(f"未找到PMID: {missing_pmid} ({missing_pmid/total_rows*100:.1f}%)")
    
    # 保存结果
    df.to_excel(output_file, index=False)
    print(f"\n结果已保存到: {output_file}")
    
    # 显示前几行预览
    print("\n前10行预览:")
    preview_df = df[[doi_column, 'PMID']].head(10)
    for idx, row in preview_df.iterrows():
        doi_short = row[doi_column].split('doi.org/')[-1] if pd.notna(row[doi_column]) and 'doi.org/' in str(row[doi_column]) else row[doi_column]
        print(f"  {doi_short} -> {row['PMID']}")
    
    return df


# 主程序
if __name__ == "__main__":
    # 设置参数 - 请修改以下三个参数
    input_file = "../../data/RetroBioCat/trial_activity_data.xlsx"  # 输入文件路径
    output_file = "../../data/RetroBioCat/trial_activity_data_with_pmid.xlsx"  # 输出文件路径
    your_email = "your_email@example.com"  # 请替换为你的实际邮箱
    
    # 处理文件
    df_result = process_excel_with_doi(
        input_file=input_file,
        output_file=output_file,
        doi_column='html_doi',
        email=your_email
    )
    
    print("\n完成！")
