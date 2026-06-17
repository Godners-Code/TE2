import os
import subprocess
import shutil
import re
from urllib.parse import urlparse

def modify_links(html_content):
    """
    使用极其严格且全面的正则，替换 HTML 标签中所有指向 .md 的超链接，
    同时确保不会误伤外部网页链接（如 http://.../file.md）或邮箱。
    """
    # 1. 核心正则 1：精准匹配 href="...path/file.md" 或 href='...path/file.md'
    # 允许链接后面带有 #锚点 或 ?参数
    pattern = r'href=(["\'])([^"\']+?)\.md(?:([?#][^"\']*))?\1'
    
    def replacer(match):
        quote = match.group(1)       # 引号 (" 或 ')
        path = match.group(2)        # 路径主干 (filename)
        suffix = match.group(3) or '' # 锚点或参数 (如 #header)，没有则为空字符串
        
        # 排除外部绝对链接（例如 http://github.com/xxx/file.md 不应该被修改）
        parsed = urlparse(path)
        if parsed.scheme or parsed.netloc:
            return match.group(0) # 原样返回，不作修改
            
        return f'href={quote}{path}.html{suffix}{quote}'

    updated = re.sub(pattern, replacer, html_content, flags=re.IGNORECASE)
    return updated

def main():
    # 获取当前工作目录的绝对路径
    base_dir = os.path.abspath(".")
    site_dir = os.path.join(base_dir, "_site")

    # 1. 初始化 _site 目录
    if os.path.exists(site_dir):
        shutil.rmtree(site_dir)
    os.makedirs(site_dir, exist_ok=True)

    # 2. 复制 Images 资源
    images_src = os.path.join(base_dir, "Images")
    if os.path.exists(images_src):
        shutil.copytree(images_src, os.path.join(site_dir, "Images"), dirs_exist_ok=True)
        print("[INFO] Images copied successfully.")

    # 3. 遍历并渲染所有 .md 文件
    print("[INFO] START CONVERSION")
    for root, dirs, files in os.walk(base_dir):
        # 排除隐藏目录和输出目录
        if "_site" in root or ".git" in root or ".github" in root:
            continue
        
        for file in files:
            if file.endswith(".md"):
                md_path = os.path.join(root, file)
                
                # 计算相对于项目根目录的相对路径
                rel_path = os.path.relpath(md_path, start=base_dir)
                rel_no_ext = os.path.splitext(rel_path)[0]
                
                # 转换为绝对的目标 HTML 路径
                dest_html = os.path.join(site_dir, rel_no_ext + ".html")
                
                # 确保目标的父目录被成功创建
                dest_parent_dir = os.path.dirname(dest_html)
                os.makedirs(dest_parent_dir, exist_ok=True)
                
                print(f"[INFO] Processing: {rel_path} -> {os.path.relpath(dest_html, base_dir)}")
                
                # 调用 pandoc 渲染
                subprocess.run(["pandoc", md_path, "-s", "--metadata", "title=Voluptus Research", "-o", dest_html], check=True)
                
                # 读取生成的 HTML 并精准替换内部链接
                with open(dest_html, "r", encoding="utf-8") as f:
                    content = f.read()
                
                # 修正后的正则表达式，精准替换 href 中的 .md 后缀
                updated_content = modify_links(content)
                
                with open(dest_html, "w", encoding="utf-8") as f:
                    f.write(updated_content)
                
                # 首页主路由兜底逻辑
                if rel_path.lower() in ["readme.md", "index.md"]:
                    shutil.copyfile(dest_html, os.path.join(site_dir, "index.html"))
                    print(f"[INFO] Set {rel_path} as site root (index.html)")
                    
    print("[INFO] CONVERSION END")

    # 4. 写入 .nojekyll 文件
    with open(os.path.join(site_dir, ".nojekyll"), "w") as f:
        f.write("")

if __name__ == "__main__":
    main()
