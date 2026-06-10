"""
序列比对算法模块
包含Needleman-Wunsch全局比对和Smith-Waterman局部比对算法
"""

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class AlignmentResult:
    """
    比对结果数据类
    """
    seq1_aligned: str
    seq2_aligned: str
    score: float
    matches: int
    mismatches: int
    gaps: int
    aligned_length: int
    similarity: float
    start_pos1: int = 0
    start_pos2: int = 0
    
    @property
    def identity(self):
        """同一性百分比"""
        if self.aligned_length == 0:
            return 0.0
        return (self.matches / self.aligned_length) * 100


def _create_score_matrix(rows: int, cols: int, gap_penalty: float, is_local: bool = False):
    """
    创建得分矩阵
    """
    matrix = [[0.0] * (cols + 1) for _ in range(rows + 1)]
    
    if not is_local:
        for i in range(rows + 1):
            matrix[i][0] = i * gap_penalty
        for j in range(cols + 1):
            matrix[0][j] = j * gap_penalty
    
    return matrix


def _traceback_global(matrix, seq1: str, seq2: str, scoring, gap_penalty: float) -> Tuple[str, str]:
    """
    全局比对回溯
    """
    i, j = len(seq1), len(seq2)
    aligned1, aligned2 = [], []
    
    while i > 0 and j > 0:
        score = matrix[i][j]
        diag = matrix[i-1][j-1] + scoring.get_score(seq1[i-1], seq2[j-1])
        up = matrix[i-1][j] + gap_penalty
        left = matrix[i][j-1] + gap_penalty
        
        if score == diag:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score == up:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    while i > 0:
        aligned1.append(seq1[i-1])
        aligned2.append('-')
        i -= 1
    while j > 0:
        aligned1.append('-')
        aligned2.append(seq2[j-1])
        j -= 1
    
    return ''.join(reversed(aligned1)), ''.join(reversed(aligned2))


def _traceback_local(matrix, seq1: str, seq2: str, scoring, gap_penalty: float, max_pos: Tuple[int, int]) -> Tuple[str, str, int, int]:
    """
    局部比对回溯
    """
    i, j = max_pos
    aligned1, aligned2 = [], []
    
    while i > 0 and j > 0 and matrix[i][j] > 0:
        score = matrix[i][j]
        diag = matrix[i-1][j-1] + scoring.get_score(seq1[i-1], seq2[j-1])
        up = matrix[i-1][j] + gap_penalty
        left = matrix[i][j-1] + gap_penalty
        
        if score == diag:
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score == up:
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    return ''.join(reversed(aligned1)), ''.join(reversed(aligned2)), i, j


def _calculate_stats(aligned1: str, aligned2: str, score: float) -> AlignmentResult:
    """
    计算比对统计信息
    """
    matches = 0
    mismatches = 0
    gaps = 0
    aligned_length = len(aligned1)
    
    for a, b in zip(aligned1, aligned2):
        if a == '-' or b == '-':
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1
    
    similarity = (matches / aligned_length) * 100 if aligned_length > 0 else 0
    
    return AlignmentResult(
        seq1_aligned=aligned1,
        seq2_aligned=aligned2,
        score=score,
        matches=matches,
        mismatches=mismatches,
        gaps=gaps,
        aligned_length=aligned_length,
        similarity=similarity
    )


def needleman_wunsch(seq1: str, seq2: str, scoring) -> AlignmentResult:
    """
    Needleman-Wunsch全局比对算法
    
    Args:
        seq1: 第一条序列
        seq2: 第二条序列
        scoring: ScoringMatrix对象
        
    Returns:
        AlignmentResult对象
    """
    m, n = len(seq1), len(seq2)
    gap_penalty = scoring.gap
    
    matrix = _create_score_matrix(m, n, gap_penalty, is_local=False)
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = scoring.get_score(seq1[i-1], seq2[j-1])
            diag = matrix[i-1][j-1] + match_score
            up = matrix[i-1][j] + gap_penalty
            left = matrix[i][j-1] + gap_penalty
            matrix[i][j] = max(diag, up, left)
    
    aligned1, aligned2 = _traceback_global(matrix, seq1, seq2, scoring, gap_penalty)
    score = matrix[m][n]
    
    return _calculate_stats(aligned1, aligned2, score)


def smith_waterman(seq1: str, seq2: str, scoring) -> AlignmentResult:
    """
    Smith-Waterman局部比对算法
    
    Args:
        seq1: 第一条序列
        seq2: 第二条序列
        scoring: ScoringMatrix对象
        
    Returns:
        AlignmentResult对象
    """
    m, n = len(seq1), len(seq2)
    gap_penalty = scoring.gap
    
    matrix = _create_score_matrix(m, n, gap_penalty, is_local=True)
    
    max_score = 0.0
    max_pos = (0, 0)
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = scoring.get_score(seq1[i-1], seq2[j-1])
            diag = matrix[i-1][j-1] + match_score
            up = matrix[i-1][j] + gap_penalty
            left = matrix[i][j-1] + gap_penalty
            matrix[i][j] = max(0, diag, up, left)
            
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)
    
    aligned1, aligned2, start1, start2 = _traceback_local(matrix, seq1, seq2, scoring, gap_penalty, max_pos)
    
    result = _calculate_stats(aligned1, aligned2, max_score)
    result.start_pos1 = start1
    result.start_pos2 = start2
    
    return result
