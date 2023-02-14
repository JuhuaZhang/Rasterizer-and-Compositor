# Rasterizer and Compositor

## Requirements
- g++
- [LodePNG](https://lodev.org/lodepng/), which already included in the file
- [ImageMagick](https://imagemagick.org/script/download.php), for image comparison

## Usage
1. Generate image
    ```sh
    make build
    make run file=$filename # no space between "=", mp1xxx.txt
    ```

2. A shell script generates and compares with the reference image
    ```sh
    ./compare.sh mp1xxx.txt
    ```

## Todo List
- [ ] too dark in rgba, maybe with srgb staff.

