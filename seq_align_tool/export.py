"""
结果导出模块
支持导出为文本格式和ALN格式
"""

import os
from .alignment import AlignmentResult
from .visualization import color_alignment


def _ensure_directory(output_path):
    """确保输出目录存在"""
    directory = os.path.dirname(output_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)


def export_alignment(result: AlignmentResult, output_path: str, 
                     fmt: str = 'txt', name1: str = 'Sequence1', 
                     name2: str = 'Sequence2', line_width: int = 60) -> None:
    """
    导出比对结果到文件
    
    Args:
        result: AlignmentResult对象
        output_path: 输出文件路径
        fmt: 导出格式 ('txt' 或 'aln')
        name1: 第一条序列名称
        name2: 第二条序列名称
        line_width: 每行序列宽度
    """
    if fmt == 'aln':
        _export_aln(result, output_path, name1, name2, line_width)
    else:
        _export_txt(result, output_path, name1, name2, line_width)


def _export_txt(result: AlignmentResult, output_path: str, 
                name1: str, name2: str, line_width: int) -> None:
    """
    导出为纯文本格式
    """
    _ensure_directory(output_path)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("序列比对结果\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"序列1: {name1}\n")
        f.write(f"序列2: {name2}\n\n")
        
        f.write(f"比对得分: {result.score:.1f}\n")
        f.write(f"相似度: {result.similarity:.2f}%\n")
        f.write(f"匹配数: {result.matches} / {result.aligned_length}\n")
        f.write(f"错配数: {result.mismatches}\n")
        f.write(f"Gap数: {result.gaps}\n\n")
        
        f.write("比对详情:\n")
        f.write("-" * 70 + "\n\n")
        
        aligned1 = result.seq1_aligned
        aligned2 = result.seq2_aligned
        
        match_line = []
        for a, b in zip(aligned1, aligned2):
            if a == '-' or b == '-':
                match_line.append(' ')
            elif a == b:
                match_line.append('|')
            else:
                match_line.append('*')
        
        match_str = ''.join(match_line)
        
        max_name_len = max(len(name1), len(name2))
        
        for i in range(0, len(aligned1), line_width):
            end = min(i + line_width, len(aligned1))
            chunk1 = aligned1[i:end]
            chunk2 = aligned2[i:end]
            chunk_match = match_str[i:end]
            
            f.write(f"{name1:<{max_name_len}}: {chunk1}\n")
            f.write(f"{'':<{max_name_len}}  {chunk_match}\n")
            f.write(f"{name2:<{max_name_len}}: {chunk2}\n\n")
        
        f.write("=" * 70 + "\n")
        f.write("图例: | = 匹配, * = 错配, 空格 = Gap\n")


def _export_aln(result: AlignmentResult, output_path: str, 
                name1: str, name2: str, line_width: int) -> None:
    """
    导出为ALN格式（CLUSTAL格式）
    """
    _ensure_directory(output_path)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("CLUSTAL O(1.2.4) multiple sequence alignment\n\n\n")
        
        aligned1 = result.seq1_aligned
        aligned2 = result.seq2_aligned
        
        match_line = []
        for a, b in zip(aligned1, aligned2):
            if a == '-' or b == '-':
                match_line.append(' ')
            elif a == b:
                match_line.append('*')
            elif is_similar_residue(a, b):
                match_line.append(':')
            else:
                match_line.append('.')
        
        match_str = ''.join(match_line)
        
        name1_clean = name1.split()[0][:30] if name1 else 'seq1'
        name2_clean = name2.split()[0][:30] if name2 else 'seq2'
        
        for i in range(0, len(aligned1), line_width):
            end = min(i + line_width, len(aligned1))
            chunk1 = aligned1[i:end]
            chunk2 = aligned2[i:end]
            chunk_match = match_str[i:end]
            
            pos1_end = i + len([c for c in chunk1 if c != '-'])
            pos2_end = i + len([c for c in chunk2 if c != '-'])
            
            f.write(f"{name1_clean:<30} {chunk1} {pos1_end}\n")
            f.write(f"{name2_clean:<30} {chunk2} {pos2_end}\n")
            f.write(f"{'':<30} {chunk_match}\n\n")


def is_similar_residue(a: str, b: str) -> bool:
    """
    判断两个氨基酸残基是否相似（保守替换）
    """
    similar_groups = [
        set('AVILM'),
        set('FYW'),
        set('KRH'),
        set('DE'),
        set('NQ'),
        set('ST'),
        set('AG'),
        set('P'),
        set('C')
    ]
    
    a = a.upper()
    b = b.upper()
    
    if a == b:
        return True
    
    for group in similar_groups:
        if a in group and b in group:
            return True
    
    return False
