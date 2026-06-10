"""
批量比对模块
支持将目录中的多个序列与参考序列进行批量比对
"""

import os
from typing import List, Dict, Callable
from .alignment import AlignmentResult


def _ensure_directory(output_path):
    """确保输出目录存在"""
    directory = os.path.dirname(output_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)


def batch_alignment(reference_seq: str, query_sequences: List[tuple],
                    align_func: Callable, scoring) -> List[Dict]:
    """
    批量比对多个序列与参考序列
    
    Args:
        reference_seq: 参考序列字符串
        query_sequences: 查询序列列表，每个元素为 (name, sequence) 元组
        align_func: 比对函数 (needleman_wunsch 或 smith_waterman)
        scoring: ScoringMatrix对象
        
    Returns:
        比对结果列表，每个元素为包含比对信息的字典
    """
    results = []
    
    for name, seq in query_sequences:
        result = align_func(reference_seq, seq, scoring)
        
        results.append({
            'name': name,
            'sequence': seq,
            'result': result,
            'score': result.score,
            'similarity': result.similarity,
            'matches': result.matches,
            'mismatches': result.mismatches,
            'gaps': result.gaps,
            'aligned_length': result.aligned_length
        })
    
    results.sort(key=lambda x: x['score'], reverse=True)
    
    return results


def generate_summary_table(results: List[Dict]) -> str:
    """
    生成批量比对的汇总表格
    
    Args:
        results: 批量比对结果列表
        
    Returns:
        格式化的汇总表格字符串
    """
    lines = []
    
    header = f"{'序号':<6}{'序列名称':<30}{'得分':<12}{'相似度(%)':<12}{'匹配数':<10}{'长度':<10}"
    lines.append(header)
    lines.append("-" * len(header))
    
    for i, result in enumerate(results, 1):
        line = (f"{i:<6}{result['name']:<30}{result['score']:<12.1f}"
                f"{result['similarity']:<12.2f}{result['matches']:<10}"
                f"{result['aligned_length']:<10}")
        lines.append(line)
    
    return '\n'.join(lines)


def export_batch_results(results: List[Dict], output_path: str, 
                         reference_name: str = "Reference",
                         algorithm: str = "global") -> None:
    """
    导出批量比对结果到文件
    
    Args:
        results: 批量比对结果列表
        output_path: 输出文件路径
        reference_name: 参考序列名称
        algorithm: 使用的算法
    """
    _ensure_directory(output_path)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(f"批量比对结果\n")
        f.write(f"参考序列: {reference_name}\n")
        f.write(f"算法: {algorithm}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(generate_summary_table(results))
        f.write("\n\n")
        f.write("=" * 80 + "\n")
        f.write("详细比对结果:\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(results, 1):
            f.write(f"\n[{i}] {result['name']}\n")
            f.write("-" * 80 + "\n")
            f.write(f"得分: {result['score']:.1f}\n")
            f.write(f"相似度: {result['similarity']:.2f}%\n")
            f.write(f"匹配: {result['matches']} | 错配: {result['mismatches']} | Gap: {result['gaps']}\n")
            f.write(f"比对长度: {result['aligned_length']}\n\n")
            
            align_result = result['result']
            aligned1 = align_result.seq1_aligned
            aligned2 = align_result.seq2_aligned
            match_line = []
            
            for a, b in zip(aligned1, aligned2):
                if a == '-' or b == '-':
                    match_line.append(' ')
                elif a == b:
                    match_line.append('|')
                else:
                    match_line.append('*')
            
            match_str = ''.join(match_line)
            
            line_width = 60
            for j in range(0, len(aligned1), line_width):
                end = min(j + line_width, len(aligned1))
                f.write(f"Ref:   {aligned1[j:end]}\n")
                f.write(f"       {match_str[j:end]}\n")
                f.write(f"Query: {aligned2[j:end]}\n\n")
