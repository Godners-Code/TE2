import os, tempfile, re, time, ssl, urllib.request
from type_alias import UrlProxy, UrlHandler, UrlReq, LInt
from runtime import log
import runtime

def _cache_dir() -> str:
    """缓存目录：系统临时目录下的 cache_dir_name。"""
    return os.path.join(tempfile.gettempdir(), runtime.CFG["site"]["cache_dir_name"])

def _normalize_proxy(proxy: str) -> str:
    """规范化代理地址；空字符串表示不使用代理。"""
    p: str = (proxy or "").strip()
    if not p: return ""
    if "://" not in p: p = "http://" + p
    return p

def _opener(proxy: str = "") -> UrlProxy:
    """构造 urllib opener；proxy 为空时显式禁用代理。"""
    ctx: ssl.SSLContext = ssl.create_default_context()
    p: str = _normalize_proxy(proxy)
    proxy_handler: UrlHandler
    if p: proxy_handler = urllib.request.ProxyHandler({"http": p, "https": p})
    else: proxy_handler = urllib.request.ProxyHandler({}) # 显式禁用环境变量中的代理
    return urllib.request.build_opener(proxy_handler, urllib.request.HTTPSHandler(context=ctx))

def fetch(name: str, timeout: int = 60, proxy: str = "") -> tuple[str|None, str]:
    """从站点抓取 <name>，缓存到系统临时目录；网络失败时回退到缓存。"""
    cache: str = _cache_dir()
    path: str = os.path.join(cache, name)
    url: str = runtime.CFG["site"]["base"] + name
    t0: float = time.time()
    try:
        req: UrlReq = urllib.request.Request(url, headers={"User-Agent": runtime.CFG["site"]["user_agent"]})
        with _opener(proxy).open(req, timeout=timeout) as r: raw = r.read()
        os.makedirs(cache, exist_ok=True)
        with open(path, "wb") as f: f.write(raw)
        log("      GET %-24s %6d bytes  %.2fs" % (name, len(raw), time.time() - t0))
        return raw.decode("utf-8", "replace"), "network"
    except Exception as e:
        log("      GET %-24s FAILED (%s: %s)" % (name, type(e).__name__, e))
    if os.path.exists(path):
        txt: str = open(path, encoding="utf-8", errors="replace").read()
        log("      CACHE %-22s %6d bytes" % (name, os.path.getsize(path)))
        return txt, "cache"
    return None, "missing"

def newest_inbound_year(readme_html) -> tuple[int, LInt]:
    """从 README 解析年份最大的 Inbound-<year>.html 链接。"""
    years: LInt = sorted({int(y) for y in re.findall(r'Inbound-(\d{4})\.html', readme_html)})
    if not years: raise SystemExit("README 中未找到任何 Inbound-<year>.html 链接")
    return years[-1], years
