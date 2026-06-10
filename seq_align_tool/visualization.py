"""
可视化输出模块
支持终端彩色高亮显示比对结果
"""

import sys


class Colors:
    """
    ANSI颜色代码
    """
    RESET = '\033[0m'
    BOLD = '\033[1m'
    
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    MAGENTA = '\033[95m'
    
    BG_GREEN = '\033[42m'
    BG_RED = '\033[41m'
    BG_YELLOW = '\033[43m'
    BG_BLUE = '\033[44m'


def supports_color():
    """
    检查终端是否支持彩色输出
    """
    if not hasattr(sys.stdout, 'isatty'):
        return False
    if not sys.stdout.isatty():
        return False
    
    if sys.platform == 'win32':
        try:
            import ctypes
            kernel32 = ctypes.windll.kernel32
            kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
            return True
        except:
            return False
    
    return True


def _colorize_segment(a, b, use_color=True):
    """
    为单个碱基对添加颜色标记
    
    Args:
        a: 序列1的碱基
        b: 序列2的碱基
        use_color: 是否使用颜色
        
    Returns:
        (colored_a, colored_b, match_char) 元组
    """
    color_supported = supports_color() and use_color
    
    if a == '-' or b == '-':
        if color_supported:
            return (f"{Colors.YELLOW}{a}{Colors.RESET}", 
                    f"{Colors.YELLOW}{b}{Colors.RESET}", 
                    ' ')
        else:
            return (a, b, ' ')
    elif a == b:
        if color_supported:
            return (f"{Colors.GREEN}{a}{Colors.RESET}", 
                    f"{Colors.GREEN}{b}{Colors.RESET}", 
                    '|')
        else:
            return (a, b, '|')
    else:
        if color_supported:
            return (f"{Colors.RED}{a}{Colors.RESET}", 
                    f"{Colors.RED}{b}{Colors.RESET}", 
                    '*')
        else:
            return (a, b, '*')


def color_alignment(result, use_color=True):
    """
    为比对结果添加颜色标记
    
    Args:
        result: AlignmentResult对象
        use_color: 是否使用颜色
        
    Returns:
        (colored_seq1, colored_seq2, colored_match) 元组
    """
    aligned1 = result.seq1_aligned
    aligned2 = result.seq2_aligned
    
    colored1 = []
    colored2 = []
    match_line = []
    
    for a, b in zip(aligned1, aligned2):
        ca, cb, mc = _colorize_segment(a, b, use_color)
        colored1.append(ca)
        colored2.append(cb)
        match_line.append(mc)
    
    return (
        ''.join(colored1),
        ''.join(colored2),
        ''.join(match_line)
    )


def _colorize_chunk(chunk1, chunk2, use_color=True):
    """
    为一段序列切片添加颜色
    
    Args:
        chunk1: 序列1的切片
        chunk2: 序列2的切片
        use_color: 是否使用颜色
        
    Returns:
        (colored_chunk1, colored_chunk2, match_chunk) 元组
    """
    colored1 = []
    colored2 = []
    match_line = []
    
    for a, b in zip(chunk1, chunk2):
        ca, cb, mc = _colorize_segment(a, b, use_color)
        colored1.append(ca)
        colored2.append(cb)
        match_line.append(mc)
    
    return (
        ''.join(colored1),
        ''.join(colored2),
        ''.join(match_line)
    )


def print_alignment(result, name1='Sequence1', name2='Sequence2', 
                    line_width=60, use_color=True, show_positions=True):
    """
    打印比对结果到终端
    
    Args:
        result: AlignmentResult对象
        name1: 第一条序列名称
        name2: 第二条序列名称
        line_width: 每行显示的宽度（碱基数量）
        use_color: 是否使用彩色
        show_positions: 是否显示位置编号
    """
    raw1 = result.seq1_aligned
    raw2 = result.seq2_aligned
    
    total_length = len(raw1)
    
    max_name_len = max(len(name1), len(name2), 8)
    
    pos1 = result.start_pos1
    pos2 = result.start_pos2
    displayed_pos1 = pos1
    displayed_pos2 = pos2
    
    for i in range(0, total_length, line_width):
        end = min(i + line_width, total_length)
        chunk_raw1 = raw1[i:end]
        chunk_raw2 = raw2[i:end]
        
        chunk_colored1, chunk_colored2, chunk_match = _colorize_chunk(
            chunk_raw1, chunk_raw2, use_color
        )
        
        if show_positions:
            start_display1 = displayed_pos1 + 1
            end_display1 = displayed_pos1 + len([c for c in chunk_raw1 if c != '-'])
            start_display2 = displayed_pos2 + 1
            end_display2 = displayed_pos2 + len([c for c in chunk_raw2 if c != '-'])
            
            displayed_pos1 += len([c for c in chunk_raw1 if c != '-'])
            displayed_pos2 += len([c for c in chunk_raw2 if c != '-'])
            
            pos_str1 = f" {start_display1:5d} - {end_display1:<5d}"
            pos_str2 = f" {start_display2:5d} - {end_display2:<5d}"
        else:
            pos_str1 = ""
            pos_str2 = ""
        
        name_pad1 = f"{name1:<{max_name_len}}: "
        name_pad2 = f"{name2:<{max_name_len}}: "
        match_pad = f"{'':<{max_name_len}}  "
        
        print(f"{name_pad1}{chunk_colored1}{pos_str1}")
        print(f"{match_pad}{chunk_match}")
        print(f"{name_pad2}{chunk_colored2}{pos_str2}")
        print()
    
    separator_len = line_width + max_name_len + 15
    print("=" * separator_len)
    if use_color and supports_color():
        print(f"{Colors.CYAN}图例说明:{Colors.RESET}")
        print(f"  {Colors.GREEN}|{Colors.RESET} = 匹配 (Match)")
        print(f"  {Colors.RED}*{Colors.RESET} = 错配 (Mismatch)")
        print(f"  {Colors.YELLOW}-{Colors.RESET} = 插入缺失 (Gap)")
    else:
        print("图例说明:")
        print("  | = 匹配 (Match)")
        print("  * = 错配 (Mismatch)")
        print("  - = 插入缺失 (Gap)")


def print_scoring_matrix(matrix, residues, title="打分矩阵"):
    """
    打印打分矩阵到终端（调试用）
    """
    print(f"\n{Colors.BOLD}{title}{Colors.RESET}")
    print("=" * (len(residues) * 5 + 10))
    
    header = "     " + "  ".join(f"{r:>3}" for r in residues)
    print(header)
    
    for i, r1 in enumerate(residues):
        row = f"{r1:>4} " + "  ".join(f"{matrix[r1][r2]:>3}" for r2 in residues)
        print(row)
    
    print()
