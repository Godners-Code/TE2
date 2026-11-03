import os
from type_alias import DataImg
import runtime

def resolve_output_path(output_dir: str, output_path: str, default_name: str) -> str:
    """解析 [output].dir 与 [output].path，兼容 Windows 与 Linux。
    - output_dir 可为相对或绝对路径，统一规范化为绝对路径；
    - output_path 为空 -> output_dir/default_name；
    - output_path 不含目录部分 -> 放入 output_dir；
    - output_path 相对且含目录 -> 相对于当前工作目录；
    - output_path 为绝对路径 -> 原样使用；
    - 无扩展名 -> 自动补默认扩展名；
    - 扩展名不在支持列表 -> 抛出 ValueError。
    """
    output_dir = os.path.abspath(os.path.expanduser(output_dir or "."))
    if output_path is None or str(output_path).strip() == "":
        return os.path.join(output_dir, default_name)

    p: str = os.path.expanduser(str(output_path).strip())
    resolved: str
    if os.path.isabs(p): resolved = os.path.normpath(p)
    elif os.path.dirname(p): resolved = os.path.abspath(p)
    else: resolved = os.path.join(output_dir, p)

    img: DataImg = runtime.CFG["image"]
    supported: tuple[str, ...] = tuple[str, ...](e.lower() for e in img["supported_exts"])
    default_ext: str = str(img["default_ext"]).lower()

    ext: str = os.path.splitext(resolved)[1].lower()
    if not ext: resolved += default_ext
    elif ext not in supported:
        raise ValueError(
            "[output].path 的扩展名 %r 不受支持；可用：%s"
            % (ext, ", ".join(supported))
        )
    return resolved
