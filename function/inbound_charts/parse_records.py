import re, html
from type_alias import DataRows, SStr, LStr

def cell_text(td: str) -> str:
    """提取表格单元格的纯文本。"""
    td = re.sub(r"(?is)<br\s*/?>", " / ", td)
    td = re.sub(r"(?s)<[^>]+>", " ", td)
    td = html.unescape(td).replace("\u00a0", " ")
    return re.sub(r"\s+", " ", td).strip()

def num(s: str) -> float | None:
    """把单元格文本转换为数字，缺失值返回 None。"""
    s = (s or "").strip()
    if s in ("", "-", "--"): return None
    return float(s) if re.match(r"^-?\d+(?:\.\d+)?$", s) else None

def parse_inbound(page_html: str) -> DataRows:
    """-> 按日期升序排列的 [{date, anastrozole_mg, letrozole_mg, ...}]"""
    rows: DataRows = []
    seen: SStr = set()
    for tr in re.findall(r"(?is)<tr>(.*?)</tr>", page_html):
        tds: LStr = re.findall(r"(?is)<td[^>]*>(.*?)</td>", tr)
        if len(tds) < 6:
            continue
        v: LStr = [cell_text(c) for c in tds]
        if not re.match(r"^\d{4}-\d{2}-\d{2}$", v[0]) or v[0] in seen:
            continue
        seen.add(v[0])
        rows.append({"date": v[0] or "",
                     "anastrozole_mg": num(v[2]) or 0,
                     "letrozole_mg": num(v[3]) or 0,
                     "testosterone_mg": (num(v[4]) if len(v) > 4 else 4) or 0,
                     "hcg_IU": (num(v[5]) if len(v) > 5 else 5) or 0})
    rows.sort(key=lambda r: r["date"])
    return rows
