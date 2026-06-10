"""
GenBank序列获取模块
支持通过GenBank accession号自动获取序列（需联网）
支持本地缓存功能，避免重复联网下载
"""

import sys
import os
import json
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET


DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser('~'), '.seq_align_cache', 'genbank')


def _get_cache_dir():
    """获取缓存目录路径"""
    cache_dir = os.environ.get('SEQ_ALIGN_CACHE_DIR', DEFAULT_CACHE_DIR)
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


def _get_cache_file(accession: str, db: str = 'nucleotide') -> str:
    """获取缓存文件路径"""
    cache_dir = _get_cache_dir()
    safe_accession = accession.replace('/', '_').replace('\\', '_')
    return os.path.join(cache_dir, f"{db}_{safe_accession}.json")


def _load_from_cache(accession: str, db: str = 'nucleotide'):
    """从缓存加载序列"""
    cache_file = _get_cache_file(accession, db)
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            return data
        except Exception as e:
            print(f"读取缓存失败: {e}", file=sys.stderr)
            return None
    return None


def _save_to_cache(accession: str, name: str, sequence: str, db: str = 'nucleotide') -> None:
    """保存序列到缓存"""
    cache_file = _get_cache_file(accession, db)
    try:
        data = {
            'accession': accession,
            'name': name,
            'sequence': sequence,
            'db': db,
            'cached_at': time.strftime('%Y-%m-%d %H:%M:%S'),
            'timestamp': time.time()
        }
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"保存缓存失败: {e}", file=sys.stderr)


def _parse_fasta_data(data: str, accession: str) -> tuple:
    """解析FASTA格式数据"""
    lines = data.strip().split('\n')
    if not lines or not lines[0].startswith('>'):
        print(f"错误: 返回的数据不是有效的FASTA格式", file=sys.stderr)
        return None, None
    
    name = lines[0][1:].strip()
    sequence = ''.join(lines[1:]).upper().replace(' ', '').replace('\t', '')
    
    if not sequence:
        print(f"错误: 序列为空", file=sys.stderr)
        return None, None
    
    return name, sequence


def _fetch_from_ncbi(accession: str, db: str = 'nucleotide') -> tuple:
    """从NCBI获取序列（不使用缓存）"""
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        params = {
            'db': db,
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text',
            'tool': 'seq_align_tool',
            'email': 'user@example.com'
        }
        
        url = f"{base_url}?{urllib.parse.urlencode(params)}"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (compatible; SeqAlignTool/1.0)'
        }
        
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                print(f"错误: HTTP状态码 {response.status}", file=sys.stderr)
                return None, None
            
            data = response.read().decode('utf-8')
            
            if not data or data.startswith('Error') or 'Error' in data[:100]:
                print(f"错误: 无法获取序列 {accession}", file=sys.stderr)
                return None, None
            
            return _parse_fasta_data(data, accession)
            
    except Exception as e:
        raise


def fetch_genbank(accession: str, use_cache: bool = True) -> tuple:
    """
    从NCBI GenBank获取核苷酸序列，支持本地缓存
    
    Args:
        accession: GenBank accession号 (如: NM_001301717, U49845)
        use_cache: 是否使用本地缓存
        
    Returns:
        (name, sequence) 元组，如果获取失败返回 (None, None)
    """
    db = 'nucleotide'
    
    if use_cache:
        cached = _load_from_cache(accession, db)
        if cached:
            print(f"从本地缓存加载: {accession}", file=sys.stderr)
            return cached['name'], cached['sequence']
    
    try:
        name, sequence = _fetch_from_ncbi(accession, db)
        if name and sequence:
            if use_cache:
                _save_to_cache(accession, name, sequence, db)
            return name, sequence
    except urllib.error.HTTPError as e:
        print(f"HTTP错误: {e.code} - {e.reason}", file=sys.stderr)
    except urllib.error.URLError as e:
        print(f"URL错误: {e.reason}", file=sys.stderr)
        print("网络连接失败，尝试使用本地缓存...", file=sys.stderr)
    except Exception as e:
        print(f"获取GenBank序列时发生错误: {e}", file=sys.stderr)
        print("尝试使用本地缓存...", file=sys.stderr)
    
    if use_cache:
        cached = _load_from_cache(accession, db)
        if cached:
            print(f"网络失败，从本地缓存加载: {accession}", file=sys.stderr)
            return cached['name'], cached['sequence']
        else:
            print("本地缓存也不存在，无法获取序列", file=sys.stderr)
    
    return None, None


def fetch_genbank_protein(accession: str, use_cache: bool = True) -> tuple:
    """
    从NCBI GenBank获取蛋白质序列，支持本地缓存
    
    Args:
        accession: GenBank蛋白质accession号
        use_cache: 是否使用本地缓存
        
    Returns:
        (name, sequence) 元组，如果获取失败返回 (None, None)
    """
    db = 'protein'
    
    if use_cache:
        cached = _load_from_cache(accession, db)
        if cached:
            print(f"从本地缓存加载: {accession}", file=sys.stderr)
            return cached['name'], cached['sequence']
    
    try:
        name, sequence = _fetch_from_ncbi(accession, db)
        if name and sequence:
            if use_cache:
                _save_to_cache(accession, name, sequence, db)
            return name, sequence
    except Exception as e:
        print(f"获取GenBank蛋白质序列时发生错误: {e}", file=sys.stderr)
        print("尝试使用本地缓存...", file=sys.stderr)
    
    if use_cache:
        cached = _load_from_cache(accession, db)
        if cached:
            print(f"网络失败，从本地缓存加载: {accession}", file=sys.stderr)
            return cached['name'], cached['sequence']
    
    return None, None


def clear_cache(accession: str = None, db: str = None) -> int:
    """
    清理缓存
    
    Args:
        accession: 要清理的accession号，如果为None则清理所有缓存
        db: 数据库类型，如果为None则清理所有数据库
        
    Returns:
        清理的文件数量
    """
    cache_dir = _get_cache_dir()
    count = 0
    
    if accession:
        if db:
            cache_file = _get_cache_file(accession, db)
            if os.path.exists(cache_file):
                os.remove(cache_file)
                count += 1
        else:
            for db_type in ['nucleotide', 'protein']:
                cache_file = _get_cache_file(accession, db_type)
                if os.path.exists(cache_file):
                    os.remove(cache_file)
                    count += 1
    else:
        if os.path.exists(cache_dir):
            for filename in os.listdir(cache_dir):
                if filename.endswith('.json'):
                    if db is None or filename.startswith(f"{db}_"):
                        os.remove(os.path.join(cache_dir, filename))
                        count += 1
    
    return count


def list_cache() -> list:
    """
    列出所有缓存的序列
    
    Returns:
        缓存条目列表，每个元素包含accession, name, db, cached_at信息
    """
    cache_dir = _get_cache_dir()
    entries = []
    
    if not os.path.exists(cache_dir):
        return entries
    
    for filename in sorted(os.listdir(cache_dir)):
        if filename.endswith('.json'):
            try:
                with open(os.path.join(cache_dir, filename), 'r', encoding='utf-8') as f:
                    data = json.load(f)
                entries.append({
                    'accession': data.get('accession'),
                    'name': data.get('name'),
                    'db': data.get('db'),
                    'cached_at': data.get('cached_at'),
                    'sequence_length': len(data.get('sequence', ''))
                })
            except:
                pass
    
    return entries


def search_genbank(term: str, db: str = 'nucleotide', max_results: int = 10) -> list:
    """
    在GenBank中搜索序列
    
    Args:
        term: 搜索关键词
        db: 数据库 ('nucleotide' 或 'protein')
        max_results: 最大返回结果数
        
    Returns:
        包含accession号和描述的列表
    """
    try:
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        
        params = {
            'db': db,
            'term': term,
            'retmax': max_results,
            'retmode': 'xml',
            'tool': 'seq_align_tool',
            'email': 'user@example.com'
        }
        
        url = f"{search_url}?{urllib.parse.urlencode(params)}"
        
        headers = {
            'User-Agent': 'Mozilla/5.0 (compatible; SeqAlignTool/1.0)'
        }
        
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                return []
            
            data = response.read().decode('utf-8')
            root = ET.fromstring(data)
            
            ids = [id_elem.text for id_elem in root.findall('.//Id')]
            
            if not ids:
                return []
            
            fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {
                'db': db,
                'id': ','.join(ids),
                'rettype': 'fasta',
                'retmode': 'text',
                'tool': 'seq_align_tool',
                'email': 'user@example.com'
            }
            
            url = f"{fetch_url}?{urllib.parse.urlencode(params)}"
            req = urllib.request.Request(url, headers=headers)
            
            with urllib.request.urlopen(req, timeout=30) as response2:
                data2 = response2.read().decode('utf-8')
                
                results = []
                entries = data2.split('>')[1:]
                
                for entry in entries[:max_results]:
                    lines = entry.split('\n')
                    if lines:
                        header = lines[0].strip()
                        accession = header.split()[0] if header else ''
                        description = header[len(accession):].strip() if accession else header
                        results.append({
                            'accession': accession,
                            'description': description
                        })
                
                return results
            
    except Exception as e:
        print(f"搜索GenBank时发生错误: {e}", file=sys.stderr)
        return []
