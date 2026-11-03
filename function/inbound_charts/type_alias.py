from datetime import datetime
from numpy.typing import NDArray
import numpy, urllib.request, argparse
type LDT    = list[datetime]
type LStr   = list[str]
type LInt   = list[int]
type LFloat = list[float]
type SStr   = set[str]
type NPFloats = NDArray[numpy.float64]
type DFloat = dict[str, float]
type DInt   = dict[str, int]
type DStr   = dict[str, str]
type DictParam = dict[str, str|int]
type DataRow   = dict[str, str|float]
type DataRows  = list[DataRow]
type DataPages = dict[int, DataRows]
type DataDose  = tuple[datetime, float|str]
type DataDoses = list[DataDose]
type DBDoses   = dict[str, DataDoses]
type DataAxes  = dict[str, float|NPFloats]
type DataImg   = dict[str, LStr|str]
type UrlProxy   = urllib.request.OpenerDirector
type UrlReq     = urllib.request.Request
type UrlHandler = urllib.request.ProxyHandler
type ArgName    = argparse.Namespace
type ArgParser  = argparse.ArgumentParser
type TConfig = dict[str, dict]
type TMeta   = dict[str, bool|str|dict]