# Rasterizer and Compositor

## Project Overview
Generates the images from the text files. 

The program can rasterizes triangles, lines and points.

- Rasterizing triangles, it use [Digital differential analyzer](https://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm))(DDA) algorithm. It supports the following features:
    - Depth buffer: Rendering objects from far to close.
    - sRGB - RGB conversion: convert colors from sRGB to linear color space before interpolating and convert back to sRGB after interpolating.
    - Frustum clipping: Perform clipping on each primitive befroe rendering.
    
 - Rasterizing lines, it use [Bresenham's line algorithm](https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm). Drawing an 8-connect line between the two given vertices.

 - Rasterizing points. Filling all pixels in a given square. Support depth buffer.

Note: The text files are in the `reference` folder.

## Requirements
- g++
- [LodePNG](https://lodev.org/lodepng/), which already included in the file

## Usage
```sh
    make build
    make run file=$filename # no space between "=", mp1xxx.txt
```

## Examples

### Basic Triangles

<img src="./generated/indexing.png"  width="5%"> <img src="./generated/tri1.png"  width="5%"> <img src="./generated/tri2.png"  width="5%"> 

### Depth Image
<img src="./generated/depth.png"  width="10%"> 

### Frustum:
<img src="./generated/frustum.png"  width="10%"> 

### sRGB
<img src="./generated/srgb.png"  width="10%"> 

### Line and Point
<img src="./generated/point.png"  width="10%"> <img src="./generated/line.png"  width="10%"> 

