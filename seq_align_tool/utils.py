"""
工具函数模块
包含DNA序列处理、反向互补等通用功能
"""

from typing import Dict, Tuple
from .alignment import AlignmentResult, needleman_wunsch, smith_waterman


COMPLEMENT_TABLE = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'N': 'N', 'n': 'n',
    'U': 'A', 'u': 'a',
    'R': 'Y', 'Y': 'R', 'r': 'y', 'y': 'r',
    'M': 'K', 'K': 'M', 'm': 'k', 'k': 'm',
    'W': 'W', 'S': 'S', 'w': 'w', 's': 's',
    'B': 'V', 'V': 'B', 'b': 'v', 'v': 'b',
    'D': 'H', 'H': 'D', 'd': 'h', 'h': 'd',
    '-': '-'
}


def reverse_complement(seq: str) -> str:
    """
    计算DNA序列的反向互补序列
    
    Args:
        seq: DNA序列字符串
        
    Returns:
        反向互补序列字符串
    """
    try:
        complement = [COMPLEMENT_TABLE.get(base, base) for base in seq]
        return ''.join(reversed(complement))
    except Exception as e:
        print(f"计算反向互补时出错: {e}")
        return seq


def align_with_reverse_complement(seq1: str, seq2: str, scoring, 
                                  algorithm: str = 'global') -> Dict:
    """
    比对DNA序列时，自动比较原序列和反向互补序列，选择得分更高的方向
    
    Args:
        seq1: 第一条DNA序列
        seq2: 第二条DNA序列
        scoring: ScoringMatrix对象
        algorithm: 比对算法 ('global' 或 'local')
        
    Returns:
        包含比对结果和方向信息的字典:
        {
            'result': AlignmentResult对象,
            'direction': 'forward' 或 'reverse_complement',
            'seq2_used': 实际使用的seq2序列（原序列或反向互补）,
            'forward_score': 原方向得分,
            'reverse_score': 反向互补得分
        }
    """
    align_func = needleman_wunsch if algorithm == 'global' else smith_waterman
    
    rc_seq2 = reverse_complement(seq2)
    
    forward_result = align_func(seq1, seq2, scoring)
    reverse_result = align_func(seq1, rc_seq2, scoring)
    
    if reverse_result.score > forward_result.score:
        return {
            'result': reverse_result,
            'direction': 'reverse_complement',
            'seq2_used': rc_seq2,
            'forward_score': forward_result.score,
            'reverse_score': reverse_result.score
        }
    else:
        return {
            'result': forward_result,
            'direction': 'forward',
            'seq2_used': seq2,
            'forward_score': forward_result.score,
            'reverse_score': reverse_result.score
        }


def is_dna_sequence(seq: str) -> bool:
    """
    检测序列是否为DNA序列
    
    Args:
        seq: 序列字符串
        
    Returns:
        True如果是DNA序列，False如果是蛋白质序列
    """
    dna_chars = set('ATCGNURYKMSWBDHV-')
    seq_upper = seq.upper().replace('-', '').replace(' ', '')
    
    if not seq_upper:
        return False
    
    protein_chars = set('ACDEFGHIKLMNPQRSTVWY*')
    dna_only = set('TUX')
    
    if any(c in dna_only for c in seq_upper):
        return True
    
    if all(c in dna_chars for c in seq_upper):
        return True
    
    if all(c in protein_chars for c in seq_upper):
        return False
    
    dna_count = sum(1 for c in seq_upper if c in 'ATCG')
    protein_count = sum(1 for c in seq_upper if c in 'EFILPQ')
    
    return dna_count > protein_count


def format_sequence(seq: str, line_width: int = 80) -> str:
    """
    格式化序列为指定宽度的多行格式
    
    Args:
        seq: 序列字符串
        line_width: 每行宽度
        
    Returns:
        格式化的序列字符串
    """
    return '\n'.join([seq[i:i+line_width] for i in range(0, len(seq), line_width)])


def calculate_gc_content(seq: str) -> float:
    """
    计算DNA序列的GC含量百分比
    
    Args:
        seq: DNA序列字符串
        
    Returns:
        GC含量百分比
    """
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    total = len([c for c in seq_upper if c in 'ATCG'])
    
    if total == 0:
        return 0.0
    
    return (gc_count / total) * 100
