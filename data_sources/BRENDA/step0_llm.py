#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step0_header_standardizer.py - 通用Header标准化工具
==================================================
使用LLM将任意CSV列名映射到Schema标准字段

适用于所有数据源: SABIO-RK, Rhea, BRENDA
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

# 添加项目根目录到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from enzymeml_utils import SchemaLoader

try:
    import yaml
except ImportError:
    yaml = None

try:
    from openai import OpenAI
    _HAS_OPENAI_SDK = True
except ImportError:
    OpenAI = None
    _HAS_OPENAI_SDK = False


# ============================================================
# I/O工具
# ============================================================

def _try_read_csv(path: Path, encoding: Optional[str] = None, sep: Optional[str] = None) -> pd.DataFrame:
    """尝试多种编码和分隔符读取CSV"""
    encs = [encoding] if encoding else ["utf-8", "utf-8-sig", "gb18030", "cp1252", "latin1"]
    seps = [sep] if sep else [None, ",", "\t", ";", "|"]
    last_err = None
    
    for enc in encs:
        for s in seps:
            try:
                return pd.read_csv(path, encoding=enc, sep=s)
            except Exception as e:
                last_err = e
                continue
    
    raise RuntimeError(f"Failed to read CSV: {path} (last error: {last_err})")


def read_table_any(path: Path, encoding: Optional[str] = None, sep: Optional[str] = None) -> pd.DataFrame:
    """读取各种格式的表格文件"""
    ext = path.suffix.lower()
    
    if ext in [".csv", ".tsv"]:
        return _try_read_csv(path, encoding=encoding, sep=sep)
    elif ext in [".parquet", ".pq"]:
        return pd.read_parquet(path)
    elif ext in [".json"]:
        df = pd.read_json(path, orient="records", lines=False)
        if not isinstance(df, pd.DataFrame):
            raise RuntimeError("JSON did not parse to a DataFrame.")
        return df
    elif ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    else:
        return _try_read_csv(path, encoding=encoding, sep=sep)


def default_output_path(input_path: Path, explicit: Optional[str] = None) -> Path:
    """生成默认输出路径"""
    if explicit:
        return Path(explicit)
    return input_path.with_name(input_path.stem + "_standardized.csv")


# ============================================================
# LLM Header映射
# ============================================================

class HeaderMapper:
    """基于LLM的Header映射器"""
    
    def __init__(
        self, 
        api_key: Optional[str],
        model: str = "gpt-4o",
        temperature: float = 0.0,
        max_tokens: int = 2500
    ):
        self.api_key = api_key
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.has_llm = api_key and _HAS_OPENAI_SDK
    
    def map_headers(
        self,
        source_headers: List[str],
        target_fields: List[str],
        verbose: bool = False
    ) -> Dict[str, Any]:
        """使用LLM将源header映射到目标schema字段"""
        if not self.has_llm:
            if verbose:
                reason = "no API key" if not self.api_key else "OpenAI SDK not installed"
                print(f"[step0] LLM disabled ({reason}); using identity mapping.")
            
            mappings = [{"source": h, "canonical": h} for h in source_headers]
            return {"mappings": mappings, "notes": "LLM disabled; identity mapping used."}
        
        prompt = self._build_prompt(source_headers, target_fields)
        
        try:
            response = self._call_llm(prompt)
            payload = self._parse_response(response)
            
            if verbose:
                mapped_count = sum(
                    1 for m in payload.get('mappings', [])
                    if m.get('canonical', '').lower() != 'unmapped'
                )
                print(f"[step0] LLM mapped {mapped_count}/{len(source_headers)} headers")
            
            return payload
            
        except Exception as e:
            if verbose:
                print(f"[step0] LLM mapping failed: {e}")
                print(f"[step0] Falling back to identity mapping")
            
            mappings = [{"source": h, "canonical": h} for h in source_headers]
            return {"mappings": mappings, "notes": f"LLM error: {e}"}
    
    def _build_prompt(self, source_headers: List[str], target_fields: List[str]) -> str:
        """构建LLM提示词"""
        # 格式化目标schema字段
        if target_fields:
            by_class: Dict[str, List[str]] = {}
            for field in target_fields:
                if '.' in field:
                    cls, attr = field.split('.', 1)
                    by_class.setdefault(cls, []).append(attr)
                else:
                    by_class.setdefault('_toplevel', []).append(field)
            
            schema_str = "TARGET SCHEMA FIELDS:\n"
            for cls in sorted(by_class.keys()):
                schema_str += f"\n{cls}:\n"
                for attr in sorted(by_class[cls]):
                    full_path = f"{cls}.{attr}" if cls != '_toplevel' else attr
                    schema_str += f"  - {full_path}\n"
        else:
            schema_str = "No schema provided. Use best judgment for field naming."
        
        source_str = "SOURCE HEADERS TO MAP:\n" + "\n".join(f"- {h}" for h in source_headers)
        
        user_prompt = f"""Map these source headers to target schema fields based on SEMANTIC MEANING.

RULES:
1. Match by MEANING, not just string similarity
   - "temp" could map to "Measurement.temperature" 
   - "enz_name" could map to "Protein.name"
   - "substrate_conc" could map to "MeasurementData.initial"

2. Use "unmapped" when:
   - The header is an internal ID with no schema equivalent
   - The meaning is unclear or ambiguous
   - No reasonable semantic match exists

3. Handle duplicates intelligently:
   - If one field clearly has multiple meanings (e.g., "modifier_1", "modifier_2"), map to same target
   - The system will add .1, .2 suffixes automatically later
   - But AVOID mapping truly different concepts to the same target

4. Output format:
   - ONLY output valid JSON
   - No explanations, no markdown formatting
   - Exact structure shown below

{schema_str}

{source_str}

Output JSON:
{{
  "mappings": [
    {{"source": "header1", "canonical": "Class.field"}},
    {{"source": "header2", "canonical": "unmapped"}},
    ...
  ]
}}
"""
        
        return user_prompt
    
    def _call_llm(self, prompt: str) -> str:
        """调用OpenAI API"""
        client = OpenAI(api_key=self.api_key)
        
        system_msg = """You are an expert at mapping messy scientific dataset headers to standardized schema fields.

Your task: Map each source header to the most semantically appropriate target field from the schema, or "unmapped" if no good match exists."""
        
        resp = client.chat.completions.create(
            model=self.model,
            temperature=self.temperature,
            max_tokens=self.max_tokens,
            messages=[
                {"role": "system", "content": system_msg},
                {"role": "user", "content": prompt},
            ],
        )
        
        return resp.choices[0].message.content or "{}"
    
    def _parse_response(self, content: str) -> Dict[str, Any]:
        """解析LLM JSON响应"""
        import re
        
        # 移除markdown包装
        if '```json' in content:
            content = content.split('```json')[1].split('```')[0]
        elif '```' in content:
            content = content.split('```')[1].split('```')[0]
        
        content = content.strip()
        
        # 提取JSON
        match = re.search(r'\{[\s\S]*\}', content)
        if match:
            return json.loads(match.group(0))
        else:
            return json.loads(content)


# ============================================================
# 映射应用
# ============================================================

class MappingApplier:
    """应用header映射，自动处理重复"""
    
    def apply(
        self,
        df: pd.DataFrame,
        mapping: Dict[str, Any],
        verbose: bool = False
    ) -> pd.DataFrame:
        """应用header映射，重复列自动添加.1, .2后缀"""
        if not mapping or "mappings" not in mapping:
            if verbose:
                print("[step0] No mapping provided; returning unchanged.")
            return df.copy()
        
        df_out = df.copy()
        rename_plan: Dict[str, str] = {}
        usage_count: Dict[str, int] = {}
        
        # 构建重命名计划
        for item in mapping.get("mappings", []):
            src_col = str(item.get("source", ""))
            target_field = item.get("canonical")
            
            if not src_col or src_col not in df_out.columns:
                continue
            
            if not target_field or str(target_field).lower() == "unmapped":
                continue
            
            new_name = self._get_unique_name(
                str(target_field),
                usage_count,
                df_out.columns,
                rename_plan
            )
            rename_plan[src_col] = new_name
        
        # 应用重命名
        if rename_plan:
            df_out = df_out.rename(columns=rename_plan)
            
            if verbose:
                print(f"[step0] Renamed {len(rename_plan)} columns:")
                for old, new in sorted(rename_plan.items())[:20]:
                    print(f"  {old:<30} -> {new}")
                if len(rename_plan) > 20:
                    print(f"  ... and {len(rename_plan) - 20} more")
        else:
            if verbose:
                print("[step0] No columns renamed")
        
        return df_out
    
    def _get_unique_name(
        self,
        base: str,
        usage_count: Dict[str, int],
        existing_cols: pd.Index,
        rename_plan: Dict[str, str]
    ) -> str:
        """生成唯一列名，如需要添加后缀"""
        if base not in usage_count and base not in existing_cols and base not in rename_plan.values():
            usage_count[base] = 0
            return base
        
        idx = usage_count.get(base, 0) + 1
        while f"{base}.{idx}" in existing_cols or f"{base}.{idx}" in rename_plan.values():
            idx += 1
        
        usage_count[base] = idx
        return f"{base}.{idx}"


# ============================================================
# 主管道
# ============================================================

class HeaderStandardizer:
    """Header标准化主协调器"""
    
    def __init__(
        self,
        schema_path: Optional[Path] = None,
        api_key: Optional[str] = None,
        model: str = "gpt-4o",
        temperature: float = 0.0,
        max_tokens: int = 2500
    ):
        self.schema_registry = SchemaLoader(schema_path)
        self.mapper = HeaderMapper(api_key, model, temperature, max_tokens)
        self.applier = MappingApplier()
    
    def standardize(
        self,
        input_path: Path,
        output_path: Path,
        encoding: Optional[str] = None,
        sep: Optional[str] = None,
        verbose: bool = False
    ):
        """运行完整的header标准化管道"""
        # 加载数据
        if verbose:
            print(f"[step0] Loading: {input_path}")
        
        df = read_table_any(input_path, encoding=encoding, sep=sep)
        source_headers = [str(c) for c in df.columns]
        
        if verbose:
            print(f"[step0] Loaded {len(df)} rows x {len(source_headers)} columns")
        
        # 获取schema中的目标字段
        target_fields = self.schema_registry.get_all_field_paths()
        
        if verbose and target_fields:
            print(f"[step0] Schema has {len(target_fields)} target fields")
        
        # 映射headers
        if verbose:
            print(f"[step0] Mapping headers via LLM...")
        
        mapping = self.mapper.map_headers(
            source_headers=source_headers,
            target_fields=target_fields,
            verbose=verbose
        )
        
        # 应用映射
        df_standardized = self.applier.apply(df, mapping, verbose=verbose)
        
        # 保存输出
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_standardized.to_csv(output_path, index=False, encoding="utf-8")
        
        # 总结
        print(f"\n✅ Success!")
        print(f"   Input:  {input_path} ({len(df)} rows x {len(source_headers)} cols)")
        print(f"   Output: {output_path} ({len(df_standardized)} rows x {len(df_standardized.columns)} cols)")
        
        if verbose:
            mapped = sum(
                1 for m in mapping.get('mappings', [])
                if m.get('canonical', '').lower() != 'unmapped'
            )
            unmapped = len(source_headers) - mapped
            print(f"\n[step0] Mapping summary:")
            print(f"   Mapped:   {mapped}")
            print(f"   Unmapped: {unmapped}")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Universal header standardization for any dataset + schema",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # With schema (recommended)
  python step0_header_standardizer.py --input data.csv --schema schemas/enzymeml-v2-extended.yaml --api-key sk-... -v
  
  # Without schema (LLM uses general knowledge)
  python step0_header_standardizer.py --input data.csv --api-key sk-... -v
  
  # Custom output
  python step0_header_standardizer.py --input data.csv --schema schemas/enzymeml-v2-extended.yaml --output clean.csv --api-key sk-...
        """
    )
    
    parser.add_argument("--input", required=True,
                       help="Input file (.csv/.xlsx/.parquet/.json/.tsv)")
    parser.add_argument("--schema",
                       help="LinkML YAML schema (optional but recommended)")
    parser.add_argument("--output",
                       help="Output path (default: <input>_standardized.csv)")
    parser.add_argument("--api-key",
                       help="OpenAI API key (or set OPENAI_API_KEY env var)")
    parser.add_argument("--model", default="gpt-4o",
                       help="OpenAI model (default: gpt-4o)")
    parser.add_argument("--temperature", type=float, default=0.0,
                       help="LLM temperature (default: 0.0)")
    parser.add_argument("--max-tokens", type=int, default=2500,
                       help="Max LLM response tokens (default: 2500)")
    parser.add_argument("--encoding",
                       help="Force input encoding (auto-detect if not set)")
    parser.add_argument("--sep",
                       help="Force CSV separator (auto-detect if not set)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # 验证输入
    input_path = Path(args.input)
    if not input_path.exists():
        raise SystemExit(f"ERROR: Input file not found: {input_path}")
    
    # 获取路径
    schema_path = Path(args.schema) if args.schema else None
    output_path = default_output_path(input_path, args.output)
    
    # 获取API key
    api_key = args.api_key or os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("WARNING: No API key provided; will use identity mapping (no changes)")
    
    # 运行标准化
    standardizer = HeaderStandardizer(
        schema_path=schema_path,
        api_key=api_key,
        model=args.model,
        temperature=args.temperature,
        max_tokens=args.max_tokens
    )
    
    standardizer.standardize(
        input_path=input_path,
        output_path=output_path,
        encoding=args.encoding,
        sep=args.sep,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()