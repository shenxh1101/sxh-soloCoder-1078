"""
保守性分析模块
支持多序列保守性分析和ASCII保守度曲线图
"""

from typing import List, Dict, Tuple
from collections import Counter
from .alignment import needleman_wunsch
from .scoring import ScoringMatrix


def _align_multiple_sequences(sequences: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """
    简单的多序列比对（使用渐进式比对策略）
    先将所有序列与第一个序列进行两两比对，然后构建多序列比对
    """
    if len(sequences) < 2:
        return sequences
    
    aligned_sequences = []
    master_name, master_seq = sequences[0]
    aligned_sequences.append((master_name, master_seq))
    
    scoring = ScoringMatrix(seq_type='dna' if _detect_seq_type(master_seq) == 'dna' else 'protein',
                           match=2, mismatch=-1, gap=-2)
    
    for name, seq in sequences[1:]:
        result = needleman_wunsch(master_seq, seq, scoring)
        aligned_sequences.append((name, result.seq2_aligned))
    
    return aligned_sequences


def _detect_seq_type(seq: str) -> str:
    """检测序列类型"""
    dna_chars = set('ATCGN')
    seq_upper = seq.upper().replace('-', '')
    if all(c in dna_chars for c in seq_upper):
        return 'dna'
    return 'protein'


def calculate_conservation(aligned_sequences: List[Tuple[str, str]]) -> Dict:
    """
    计算多序列比对每个位置的保守度
    
    Args:
        aligned_sequences: 比对后的序列列表，每个元素为 (name, aligned_seq) 元组
        
    Returns:
        包含保守性信息的字典
    """
    if len(aligned_sequences) < 2:
        raise ValueError("保守性分析需要至少2个序列")
    
    align_length = len(aligned_sequences[0][1])
    num_sequences = len(aligned_sequences)
    
    conservation_scores = []
    conservation_levels = []
    
    for pos in range(align_length):
        residues = []
        for _, seq in aligned_sequences:
            if pos < len(seq):
                residues.append(seq[pos])
        
        count = Counter(residues)
        
        if '-' in count:
            gap_count = count['-']
            if gap_count == num_sequences:
                conservation_scores.append(0.0)
                conservation_levels.append('gap')
                continue
        
        non_gap_residues = [r for r in residues if r != '-']
        
        if not non_gap_residues:
            conservation_scores.append(0.0)
            conservation_levels.append('gap')
            continue
        
        most_common, most_common_count = count.most_common(1)[0]
        
        if most_common == '-':
            most_common, most_common_count = count.most_common(2)[1] if len(count) > 1 else (None, 0)
        
        conservation = (most_common_count / len(non_gap_residues)) * 100 if len(non_gap_residues) > 0 else 0
        conservation_scores.append(conservation)
        
        if conservation == 100:
            conservation_levels.append('conserved')
        elif conservation >= 80:
            conservation_levels.append('high')
        elif conservation >= 50:
            conservation_levels.append('moderate')
        else:
            conservation_levels.append('low')
    
    avg_conservation = sum(conservation_scores) / len(conservation_scores) if conservation_scores else 0
    conserved_sites = sum(1 for s in conservation_scores if s == 100)
    
    return {
        'positions': list(range(1, align_length + 1)),
        'conservation_scores': conservation_scores,
        'conservation_levels': conservation_levels,
        'avg_conservation': avg_conservation,
        'conserved_sites': conserved_sites,
        'total_sites': align_length,
        'num_sequences': num_sequences
    }


def generate_ascii_plot(conservation_scores: List[float], plot_width: int = 80, 
                       height: int = 10) -> str:
    """
    生成保守度的ASCII曲线图
    
    Args:
        conservation_scores: 保守度分数列表
        plot_width: 图表宽度（字符数）
        height: 图表高度（行数）
        
    Returns:
        ASCII图表字符串
    """
    if not conservation_scores:
        return "无数据"
    
    n = len(conservation_scores)
    
    if n > plot_width:
        bin_size = n // plot_width
        binned_scores = []
        for i in range(0, n, bin_size):
            chunk = conservation_scores[i:i+bin_size]
            if chunk:
                binned_scores.append(sum(chunk) / len(chunk))
        conservation_scores = binned_scores
        n = len(conservation_scores)
    
    plot_lines = []
    
    for h in range(height, -1, -1):
        threshold = h * (100.0 / height)
        line = '|'
        
        for score in conservation_scores:
            if score >= threshold:
                if score >= 100:
                    line += '*'
                elif score >= 80:
                    line += '#'
                elif score >= 50:
                    line += '+'
                else:
                    line += '.'
            else:
                line += ' '
        
        if h % 2 == 0:
            line += f' {int(threshold):3d}%'
        plot_lines.append(line)
    
    x_axis = '+' + '-' * len(conservation_scores) + '+'
    plot_lines.append(x_axis)
    
    if n > 0:
        label_line = '|'
        step = max(1, n // 5)
        for i in range(n):
            if i % step == 0:
                label_line += str((i + 1) % 10)
            else:
                label_line += ' '
        label_line += '| Position'
        plot_lines.append(label_line)
    
    legend = []
    legend.append("")
    legend.append("图例: * = 100% (完全保守), # = >=80%, + = >=50%, . = <50%")
    legend.append(f"平均保守度: {sum(conservation_scores)/len(conservation_scores):.1f}%")
    legend.append(f"完全保守位点: {sum(1 for s in conservation_scores if s >= 100)} / {len(conservation_scores)}")
    
    plot_lines.extend(legend)
    
    return '\n'.join(plot_lines)


def conservation_analysis(sequences: List[Tuple[str, str]], plot_width: int = 80) -> Dict:
    """
    执行完整的保守性分析
    
    Args:
        sequences: 序列列表，每个元素为 (name, sequence) 元组
        plot_width: ASCII曲线图宽度
        
    Returns:
        包含保守性分析结果的字典
    """
    if len(sequences) < 2:
        raise ValueError("保守性分析需要至少2个序列")
    
    aligned = _align_multiple_sequences(sequences)
    
    result = calculate_conservation(aligned)
    
    ascii_plot = generate_ascii_plot(result['conservation_scores'], plot_width=plot_width)
    result['ascii_plot'] = ascii_plot
    result['aligned_sequences'] = aligned
    
    return result


def print_conservation_result(result: Dict) -> None:
    """
    打印保守性分析结果
    """
    print("\n" + "="*70)
    print("多序列保守性分析结果")
    print("="*70)
    print(f"\n序列数量: {result['num_sequences']}")
    print(f"比对长度: {result['total_sites']} 个位置")
    print(f"平均保守度: {result['avg_conservation']:.2f}%")
    print(f"完全保守位点: {result['conserved_sites']} / {result['total_sites']} "
          f"({result['conserved_sites']/result['total_sites']*100:.1f}%)")
    print(f"\n保守度曲线图:")
    print("-" * 70)
    print(result['ascii_plot'])
    print("=" * 70)
