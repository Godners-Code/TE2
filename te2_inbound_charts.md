<!-- markdownlint-disable MD060 -->
# 使用说明 — te2_inbound_charts

**T/E2 ADJUST** 的一站式入站处理脚本。脚本从 TE2 站点抓取每日给药记录，运行一室一级吸收 PK 模型，并绘制最终的双轴图表。

**所有可调参数均在 `te2_inbound_charts.toml` 中配置**，命令行可用 `--config` 显式指定配置文件路径）。

---

## 1. 环境要求

| 项目 | 要求 |
|:----:|:----|
| Python | **3.12 或更高版本**（低于此版本脚本会直接退出） |
| 操作系统 | Windows / Linux / macOS |
| 第三方库 | `matplotlib`、`numpy`、`scipy` |
| 网络 | 能访问站点基址（如 <https://te2.tianyue.ren/>），直连或经代理 |

安装依赖：

```bash
python -m pip install matplotlib numpy scipy
```

> tomllib 是 Python 3.11+ 内置模块，无需单独安装。

## 2. 目录结构

脚本依赖同目录下的配置文件，命名规则为 **脚本名 +** .toml：

```plaintext
te2_inbound_charts.py
te2_inbound_charts.toml
USAGE.md
```

若从其他目录运行，请使用 `--config` 显式指向 TOML 文件的路径：

```bash
python /path/to/te2_inbound_charts.py --config /path/to/custom_configure.toml
```

## 3. 快速开始

```bash
# 默认运行：最后记录日、180 天窗口、输出到当前目录，使用脚本同目录的 TOML 配置文件
python te2_inbound_charts.py
```

默认输出：

```plaintext
./Inbound-<结束日期yyyymmdd>-<统计天数>.png
./Inbound-<结束日期yyyymmdd>-<统计天数>.json
```

示例：./Inbound-20250131-180.png 以及同目录下的 ./Inbound-20250131-180.json。

元数据 JSON 始终与图片**同名、同目录**，仅扩展名替换为 .json（如 chart.jpg → chart.json，chart.webp → chart.json）。

## 4. 命令行参数

| 参数 | 默认值 | 说明 |
|:----:|:----:|:----|
| --config PATH | <脚本同目录>/<脚本名去掉扩展名>.toml | 显式指定配置文件路径 |

- 不传 --config 时，脚本自动读取同目录下的 te2_inbound_charts.toml。

## 5. 配置文件 te2_inbound_charts.toml

所有可调参数都在这里；下面按段说明。默认文件名由脚本名决定（脚本名 + .toml）。

### 5.1 [app] - 图表文字

| 键 | 说明 |
|:----:|:----|
| title | 图表标题 |
| x_label | X 轴标题 |
| ylabel_left | 左轴标题 |
| ylabel_right | 右轴标题 |
| label_anastrozole | 阿那曲唑图例文字 |
| label_letrozole | 来曲唑图例文字|

### 5.2 [site] - 站点与缓存

| 键 | 说明 |
|:----:|:----|
| base | 站点基址（如 <https://te2.tianyue.ren/>） |
| cache_dir_name | 缓存目录名；实际路径在系统临时目录下 |
| user_agent | 抓取时使用的 User-Agent 字符串 |
| proxy | 代理地址；留空表示不使用代理 |

proxy 行为：

- 为空：显式禁用代理，忽略环境变量 `HTTP_PROXY` / `HTTPS_PROXY`。  
- 非空：作为 HTTP 与 HTTPS 代理使用；未写协议前缀时会自动补 `http://`。

### 5.3 [run] - 运行窗口参数

| 键 | 默认值 | 说明 |
|:----:|:----:|:----|
| end | auto | auto = 最后记录给药日；也可填 YYYY-MM-DD |
| period | half | 统计周期 |
| window | 20 | 波动率滚动窗口天数 |

- period 取值

  | 取值 | 天数 |
  |:----:|:----|
  | monthly | 30 |
  | quarter | 90 |
  | half | 180 |
  | year | 360 |
  | <数字> | 该数字本身，如 120 |

  > 不区分大小写（Half、HALF、half 等价）。  
  > 必须满足 min_window_days <= window <= period <= max_days（默认 366）。

### 5.4 [output] — 输出位置

| 键 | 默认值 | 说明 |
|:----:|:----:|:----|
| dir | . | 输出目录；相对或绝对路径 |
| path | 空 | 输出文件路径；留空使用默认文件名 |

- path 路径解析规则：

1. 为空 → <output_dir>/Inbound-<结束日期yyyymmdd>-<统计天数>.<default_ext>
2. 绝对路径 → 原样使用。
3. 含目录部分（如 logs/te2/chart.png）→ 相对于当前工作目录，转换为绝对路径。
4. 仅文件名（如 chart.png）→ 放入 output_dir。
5. 无扩展名（如 chart）→ 自动补 .png（由 image.default_ext 决定）。
6. 扩展名不在 `image.supported_exts` 中 → 直接报错退出。

- dir 本身可为相对或绝对路径，脚本会统一转换为绝对路径。

### 5.5 [render] - 数值与采样

| 键 | 说明 |
|:----:|:----|
| sample_hour | 给药 / 采样时钟（小时） |
| sample_step_h | 日采样间隔 |
| sub_per_day | 平滑曲线相对真实采样的加密倍率 |
| n_intervals | 两个 y 轴的等分间隔数 |
| conc_decimals | 浓度规范化保留的小数位 |
| dose_lookback_h | 剂量回溯窗口（默认 960 h ≈ 40 天） |

### 5.6 [colors] - 配色

| 键 | 说明 |
|:----:|:----|
| ana_line | 阿那曲唑曲线颜色 |
| let_line | 来曲唑曲线颜色 |
| ana_bar | 阿那曲唑色带颜色 |
| let_bar | 来曲唑色带颜色 |
| grid | 网格颜色 |
| grid_lw | 网格线宽 |
| axis_text_left | 左轴文字颜色 |
| axis_text_right | 右轴文字颜色 |
| spine | 坐标轴边框颜色 |

### 5.7 [fonts] - 字体回退顺序

- monospace 为字符串数组，按顺序取第一个系统中可用的字体；
- 全部缺失时回退到 matplotlib 内置的 DejaVu Sans Mono。

### 5.8 [image] - 图片格式

| 键 | 说明 |
|:----:|:----|
| supported_exts | output_path 允许的扩展名列表 |
| default_ext | 未指定扩展名时使用的默认扩展名 |

- 支持的图片格式：

  | 扩展名 | 输出格式 |
  |:----:|:----|
  | .png（默认） | PNG |
  | .jpg、.jpeg | JPEG |
  | .bmp | BMP |
  | .webp | WebP |

- 其他扩展名（如 .gif、.svg、.tif）会收到明确的报错信息并终止运行。

### 5.9 [period] - 统计周期与边界

| 键 | 说明 |
| monthly | 天数映射（30） |
| quarter | 天数映射（90） |
| half | 天数映射（180） |
| year | 天数映射（360） |
| max_days | 统计天数上限（默认 366） |
| min_window_days | window 的最小值（默认 2） |

### 5.10 [pk.anastrozole] / [pk.letrozole] - PK 参数

| 键 | 说明 |
|:----:|:----|
| t_half_h | 半衰期（小时） |
| ka | 一级吸收速率常数（h⁻¹） |
| V_per_F | 表观分布容积（L）；设为 0 表示用标定参数反算 |
| calib_trough | 标定谷值（ng/mL） |
| calib_dose | 标定时每日剂量（mg） |
| calib_tau | 给药间隔（小时） |

- 两条路径：
  - 谷值反算（阿那曲唑采用）：V_per_F = 0，同时提供 calib_trough / calib_dose / calib_tau。
  - 直接给定（来曲唑采用）：V_per_F > 0，不需要任何 calib_*。
- 启动时校验：
  - t_half_h、ka 缺失 → 报错，配置文件
    > [pk.xxx] 缺少必需参数 t_half_h。
  - V_per_F = 0（或缺失）但缺少任一 calib_* → 报错
    > 配置文件 [pk.xxx] 的 V_per_F 为 0（表示按谷值反算），但缺少标定参数：...。
  - V_per_F > 0 时即使带上 calib_* 也不会被使用（不报错，静默忽略）。

## 6. 使用示例

```bash
# 默认运行（读取同目录的 te2_inbound_charts.toml）
python te2_inbound_charts.py

# 指定其它配置文件
python te2_inbound_charts.py --config /data/te2/te2_inbound_charts.toml

# 从其他目录运行脚本
python /path/to/te2_inbound_charts.py --config /path/to/te2_inbound_charts.toml
```

**想改变运行行为，请修改 TOML 配置后重新运行。** 例如：

- 只画 90 天窗口 → 把 [run].period 改为 `quarter`；
- 输出到固定目录 → 把 [output].dir 改为 `reports`；
- 输出为 JPEG → 把 [output].path 改为 `te2.jpg`；
- 走代理 → 把 [site].proxy 改为 `http://192.168.1.1:8080`。

以上每种情况都会在与图片同名的位置生成 .json 元数据  
（如 te2.jpg → te2.json，2025/full.png → 2025/full.json）。

## 7. 脚本做了什么

```plaintext
[1/6] 抓取 README.html -> 解析最新的 Inbound-<year>.html
[2/6] 抓取并解析每日记录（剂量、睾酮、HCG 等）
[3/6] 向前扩展窗口；必要时补入上一年度页面
[4/6] PK 仿真，每日 08:00 采样一次（给药前谷值）
[5/6] 滚动 N 天波动率 -> 平滑得到色带顶部曲线
[6/6] 绘图 + 写出与图片同名的 JSON 元数据
```

### 7.1 PK 模型

- 一室模型、一级吸收、线性叠加。
- 阿那曲唑：t½ = 50 h，ka = 2.5 h⁻¹，V/F 按稳态谷值 25.7 ng/mL（1 mg/日）标定。
- 来曲唑：t½ = 48 h，ka = 2.0 h⁻¹，V/F = 133 L（1.9 L/kg × 70 kg）。
- 剂量回溯窗口为 40 天（约 20 个半衰期），更早的剂量忽略不计。

### 7.2 图表

- **左轴**：血药浓度（ng/mL）——阿那曲唑（蓝）与来曲唑（绿），用三次样条穿过每个真实日采样点。
- **右轴**：滚动波动率（标准差，ng/mL）——堆叠色带，顶部为平滑拟合曲线，替代旧版的水平条顶段。
- **网格**：竖线在每月 1 日与 16 日；横线在两个 y 轴上间隔数相同（默认 6）。

### 7.3 字体

```plaintext
Consolas → Courier New → DejaVu Sans Mono → Liberation Mono →（matplotlib 默认等宽字体）
```

DejaVu Sans Mono 随 matplotlib 一起安装，因此即使系统中没有前几种等宽字体，图表也能正常渲染。

## 8. 缓存

```plaintext
<tempdir>/<site.cache_dir_name>/
```

- Windows：%TEMP%\te2-site\
- Linux/macOS：/tmp/te2-site/（或 $TMPDIR/te2-site/）

网络请求失败时，脚本会回退到缓存副本，并打印 `CACHE <name>` ...。如需强制重新抓取，删除该目录即可。

## 9. 输出文件

| 文件 | 说明 |
|:----:|:----|
| `<output.path>` | 生成的图片，格式由扩展名决定（默认 PNG）。 |
| `<basename>.json` | 元数据 JSON，与图片同名同目录。 |

主要元数据字段：

- file：图片文件名（例如 Inbound-20250131-180.png）。
- generated_by：生成脚本的文件名（本脚本为 te2_inbound_charts.py）。
- source.inbound_pages：使用了哪些 Inbound-`<year>`.html 页面。
- window：起止日期、天数、采样方式。
- left_axis / right_axis：标签、上限、颜色、峰值。
- volatility：滚动窗口、原始与绘制均值、样条拟合诊断。
- axis_intervals：两个坐标轴的间隔数与步长。
- pk：本次运行实际使用的 PK 参数。

## 10. 常见问题

| 现象 | 原因 / 处理 |
|:----|:----|
| Python 3.12 或更高版本是必需的。 | 升级 Python 版本。 |
| 配置文件不存在：... | 缺少 te2_inbound_charts.toml；确认脚本名与配置文件名的对应关系（脚本名去掉扩展名 + .toml），或用 --config 指向正确路径。 |
| 配置文件解析失败：... | TOML 语法错误，检查文件内容。 |
| 配置文件缺少段 [...] / 缺少 [pk.xxx] 段 | 按第 5 节补全对应段；注意需要 [run] 和 [output]，没有 [cli]。 |
| 配置文件 [pk.xxx] 缺少必需参数 t_half_h / 缺少必需参数 ka | 该药物缺少半衰期或吸收速率常数。 |
| 配置文件 [pk.xxx] 的 V_per_F 为 0（表示按谷值反算），但缺少标定参数：... | 补上 calib_trough / calib_dose / calib_tau，或直接把 V_per_F 设为正数。 |
| 配置错误：[run].period 必须为 ... | 检查 [run].period，使用 monthly / quarter / half / year 或整数天数。 |
| 配置错误：[run].period 天数必须 >= [run].window | 增大 period 或减小 window。 |
| 配置错误：[run].period 天数必须 <= <max_days> | 使用更小的统计天数，或调整 [period].max_days。 |
| 配置错误：[run].window 必须 >= <min_window_days> | 增大 window 或调整 [period].min_window_days。 |
| README.html 不可用（网络与缓存均失败） | 网络异常且无缓存；检查网络或配置 [site].proxy。 |
| Inbound-`<year>`.html 不可用 | 服务器上缺少该页面，且本地无缓存。 |
| GET ... FAILED (...) 但仍生成了图表 | 网络失败，已使用缓存；删除 `<tempdir>/te2-site/` 可强制刷新。 |
| [output].path 的扩展名 '...' 不受支持 | 改用 image.supported_exts 中的扩展名。 |
| 图表正常但中文显示为方块 | 图表按设计使用等宽字体，属预期行为。 |
| 写输出文件时提示 Permission denied | 目标目录不可写；调整 [output].dir / [output].path。 |
