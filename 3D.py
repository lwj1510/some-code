import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import xarray as xr
from metpy.units import units
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# 全局字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']


# --------------------
# 数据加载
# --------------------
def load_real_data(z_file, u_file, v_file):
    # 打开文件
    ds_z = xr.open_dataset(z_file).isel(valid_time=0)
    ds_u = xr.open_dataset(u_file).isel(valid_time=0)
    ds_v = xr.open_dataset(v_file).isel(valid_time=0)
    
    # z 转换为位势高度 (m)
    g = 9.80665
    z = ds_z['z'] / g  # 单位：m
    
    u = ds_u['u'] #.sel(valid_time=ds_u.valid_time.values[0])
    v = ds_v['v'] #.sel(valid_time=ds_v.valid_time.values[0])
    
    lons = ds_z.longitude.values
    lats = ds_z.latitude.values
    
    return lons, lats, z, u, v

# --------------------
# 三维绘图
# --------------------
def plot_circulation(lons, lats, z, u, v, level):
    fig = plt.figure(figsize=(12, 8))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1, 1, 1, projection=proj)

    # 提取该层数据
    hgt = z.sel(pressure_level=level).values
    u_wind = u.sel(pressure_level=level).values
    v_wind = v.sel(pressure_level=level).values

    # 等值线（高度场）
    cs = ax.contour(lons, lats, hgt, levels=20, colors='black', linewidths=1.0, transform=proj)
    ax.clabel(cs, inline=True, fontsize=10, fmt='%d')

    # 风场流线
    strm = ax.streamplot(lons, lats, u_wind, v_wind, density=2, color='blue', transform=proj)

    # 底图
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.set_title(f"{level} hPa 位势高度场（等值线）+ 风场（流线）", fontsize=14)

    plt.show()

# --------------------
# 主程序
# --------------------
if __name__ == "__main__":
    file_z = r"E:\Data_wrangling\data\data\ERA5\geopotential_stream-oper_daily-mean.nc"
    file_u = r"E:\Data_wrangling\data\data\ERA5\u_component_of_wind_0_daily-mean.nc"
    file_v = r"E:\Data_wrangling\data\data\ERA5\v_component_of_wind_0_daily-mean.nc"
    
    lons, lats, z, u, v = load_real_data(file_z, file_u, file_v)
    plot_circulation(lons, lats, z, u, v)
