import argparse
import csv
import re
import json
import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import geopandas as gpd
import contextily as ctx
from shapely.geometry import Point
from collections import OrderedDict
import matplotlib.colors as mcolors
from pathlib import Path
import imageio.v2 as imageio
from PIL import Image
from tqdm import tqdm

def parse_srt_file(srt_file):
    data = []
    frame_pattern = re.compile(r"FrameCnt: (\d+),")
    datetime_pattern = re.compile(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}")
    fields_pattern = re.compile(r"\[([a-zA-Z_]+): ([\d\.\-\/]+)\]")
    rel_abs_alt_pattern = re.compile(r"\[rel_alt: ([\d\.\-]+) abs_alt: ([\d\.\-]+)\]")

    with open(srt_file, 'r') as file:
        content = file.read()

    entries = content.split('\n\n')
    all_keys = set()
    for entry in entries:
        frame_info = OrderedDict()
        frame_match = frame_pattern.search(entry)
        datetime_match = datetime_pattern.search(entry)
        fields_matches = fields_pattern.findall(entry)
        rel_abs_alt_match = rel_abs_alt_pattern.search(entry)

        if frame_match:
            frame_info['FrameCnt'] = frame_match.group(1)
        if datetime_match:
            frame_info['Datetime'] = datetime_match.group(0)
        if rel_abs_alt_match:
            frame_info['rel_alt'] = rel_abs_alt_match.group(1)
            frame_info['abs_alt'] = rel_abs_alt_match.group(2)
            all_keys.add('rel_alt')
            all_keys.add('abs_alt')
        for field, value in fields_matches:
            if field not in ['rel_alt', 'abs_alt']:
                frame_info[field] = value
                all_keys.add(field)

        if frame_info:
            data.append(frame_info)

    print("Found metadata keys:", all_keys)
    return data

def write_csv(data, output_file):
    if not data:
        print("No data to write.")
        return
    
    # Define the logical order of columns
    logical_order = ['FrameCnt', 'Datetime', 'latitude', 'longitude', 'rel_alt', 'abs_alt', 'iso', 'shutter', 'fnum', 'ev', 'color_md', 'focal_len', 'ct']
    keys = logical_order + [key for key in data[0].keys() if key not in logical_order]

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        writer.writerows(data)

def write_geojson(data, output_file):
    coordinates = []
    for entry in data:
        coordinates.append([
            float(entry.get('longitude', 0)),
            float(entry.get('latitude', 0)),
            float(entry.get('abs_alt', 0))
        ])
    
    feature = {
        "type": "Feature",
        "properties": {
            "start_time": data[0].get('Datetime'),
            "end_time": data[-1].get('Datetime'),
            "start_latitude": data[0].get('latitude'),
            "start_longitude": data[0].get('longitude'),
            "end_latitude": data[-1].get('latitude'),
            "end_longitude": data[-1].get('longitude')
        },
        "geometry": {
            "type": "LineString",
            "coordinates": coordinates
        }
    }

    geojson = {
        "type": "FeatureCollection",
        "features": [feature]
    }

    with open(output_file, 'w') as geojsonfile:
        json.dump(geojson, geojsonfile, indent=4)

def create_geodataframe(data):
    geometry = [Point(float(d['longitude']), float(d['latitude']), float(d['abs_alt'])) for d in data]
    gdf = gpd.GeoDataFrame(data, geometry=geometry)
    gdf.set_crs(epsg=4326, inplace=True)  # WGS 84
    return gdf

def get_utm_crs(lon, lat):
    utm_band = str((int((lon + 180) / 6) % 60) + 1)
    if lat >= 0:
        return f"326{utm_band.zfill(2)}"
    else:
        return f"327{utm_band.zfill(2)}"

def plot_flight_summary(gdf, output_dir, srt_file, data):
    start_point = gdf.iloc[0].geometry
    end_point = gdf.iloc[-1].geometry
    start_time = data[0]['Datetime']
    end_time = data[-1]['Datetime']
    runtime = pd.to_datetime(end_time) - pd.to_datetime(start_time)
    runtime_minutes = runtime.total_seconds() / 60
    start_time_12hr = pd.to_datetime(start_time).strftime('%I:%M:%S %p')
    
    times = pd.to_datetime(gdf['Datetime'])
    norm = mcolors.Normalize(vmin=times.min().timestamp(), vmax=times.max().timestamp())
    cmap = plt.colormaps['viridis']

    # Transform to UTM for bottom plots
    utm_crs = get_utm_crs(gdf.geometry.x.mean(), gdf.geometry.y.mean())
    gdf_utm = gdf.to_crs(utm_crs)
    start_point_utm = gpd.GeoSeries([start_point], crs=gdf.crs).to_crs(utm_crs).geometry.iloc[0]
    end_point_utm = gpd.GeoSeries([end_point], crs=gdf.crs).to_crs(utm_crs).geometry.iloc[0]

    # Save top-down view
    fig, ax_top = plt.subplots(figsize=(7.5, 10))
    fig.suptitle(f"Filename: {Path(srt_file).name}\nDate: {start_time.split()[0]}\nRuntime: {runtime_minutes:.2f} min\nStart Time: {start_time_12hr}", fontsize=14)
    gdf.plot(ax=ax_top, markersize=5, color='white')
    xmin, xmax = ax_top.get_xlim()
    ymin, ymax = ax_top.get_ylim()
    ctx.add_basemap(ax_top, crs=gdf.crs.to_string(), source=ctx.providers.Esri.WorldImagery, reset_extent=False, attribution="")
    ax_top.set_xlim(xmin, xmax)
    ax_top.set_ylim(ymin, ymax)
    ax_top.scatter(start_point.x, start_point.y, color='green', s=200, label='Start', marker='x', zorder=2, edgecolor=None, linewidth=5)
    ax_top.scatter(end_point.x, end_point.y, color='red', s=200, label='End', marker='x', zorder=2, edgecolor=None, linewidth=5)
    ax_top.set_xlabel('Longitude')
    ax_top.set_ylabel('Latitude')
    ax_top.set_title('Top-Down View')
    points = gdf.geometry.apply(lambda geom: [geom.x, geom.y]).tolist()
    ax_top.scatter(*zip(*points), c=times.apply(lambda x: x.timestamp()), cmap=cmap, norm=norm, s=10, edgecolor='none', label='Flight Path')
    ax_top.legend()

    plt.tight_layout()
    top_down_file = output_dir / "top_down_view.png"
    plt.savefig(top_down_file)
    plt.close(fig)

    # Create frames for the rotating 3D plot
    frames = []
    # for angle in range(0, 360, 5):
    for angle in tqdm(range(0, 360, 5)):
        fig = plt.figure(figsize=(15, 10))
        ax_3d = fig.add_subplot(1, 1, 1, projection='3d')
        ax_3d.scatter(gdf_utm.geometry.x, gdf_utm.geometry.y, gdf_utm.geometry.z, c=times.apply(lambda x: x.timestamp()), cmap=cmap, norm=norm, s=10, edgecolor='none')
        ax_3d.scatter(start_point_utm.x, start_point_utm.y, start_point_utm.z, color='green', s=100, label='Start', marker='x', edgecolor=None, linewidth=5)
        ax_3d.scatter(end_point_utm.x, end_point_utm.y, end_point_utm.z, color='red', s=100, label='End', marker='x', edgecolor=None, linewidth=5)
        ax_3d.set_xlabel('UTM Easting (m)')
        ax_3d.set_ylabel('UTM Northing (m)')
        ax_3d.set_zlabel('Altitude (m)')
        ax_3d.set_title('3D Flight Path')
        ax_3d.view_init(elev=30, azim=angle)
        plt.tight_layout()
        
        frame_file = output_dir / f"frame_{angle}.png"
        plt.savefig(frame_file)
        plt.close(fig)

        # Merge the top-down view with the 3D plot
        top_down_img = Image.open(top_down_file)
        frame_img = Image.open(frame_file)
        combined_img = Image.new('RGB', (top_down_img.width + frame_img.width, top_down_img.height))
        combined_img.paste(top_down_img, (0, 0))
        combined_img.paste(frame_img, (top_down_img.width, 0))
        
        combined_frame_file = output_dir / f"combined_frame_{angle}.png"
        combined_img.save(combined_frame_file)
        frames.append(imageio.imread(combined_frame_file))
        os.remove(frame_file)
        os.remove(combined_frame_file)

    # Create GIF
    gif_output_file = output_dir / f"{Path(srt_file).stem}_flight_summary.gif"
    imageio.mimsave(gif_output_file, frames, fps=5)
    
    # Clean up frame files
    os.remove(top_down_file)
    
    print(f"Flight summary plot saved as {gif_output_file}")

def main():
    parser = argparse.ArgumentParser(description="Parse DJI SRT file and extract frame information.")
    parser.add_argument('srt_file', type=str, help="Path to the SRT file.")
    args = parser.parse_args()

    srt_file = Path(args.srt_file).resolve()
    base_name = srt_file.stem
    output_dir = srt_file.parent / base_name
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_output_file = output_dir / f"{base_name}_frame_index.csv"
    geojson_output_file = output_dir / f"{base_name}_flight_path.geojson"

    data = parse_srt_file(srt_file)
    write_csv(data, csv_output_file)
    write_geojson(data, geojson_output_file)
    
    gdf = create_geodataframe(data)
    plot_flight_summary(gdf, output_dir, srt_file, data)
    
    print(f"Data written to {csv_output_file} and {geojson_output_file}")

if __name__ == "__main__":
    main()
