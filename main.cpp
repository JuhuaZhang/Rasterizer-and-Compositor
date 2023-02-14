#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <math.h> // ceil()
#include <cmath>
#include "lodepng.h"

using namespace std;

int width = 0;
int height = 0;

double clip_planes[6][4] = {
    {1, 0, 0, 1},
    {-1, 0, 0, 1},
    {0, 1, 0, 1},
    {0, -1, 0, 1},
    {0, 0, 1, 1},
    {0, 0, -1, 1}};

// create vector of vertex
class vertex
{
public:
    double x;
    double y;

    double z;
    double w;

    double r;
    double g;
    double b;

    double depth;

    double a = 255; // alpha

    // empty constructor
    vertex() {}

    // xyzw + rgb, for vertex in the real world
    vertex(double _x, double _y, double _z, double _w, double _r, double _g, double _b) : x(_x), y(_y), z(_z), w(_w), r(_r), g(_g), b(_b) {}

    // xyzw + rgba, for vertex in the real world
    vertex(double _x, double _y, double _z, double _w, double _r, double _g, double _b, double _a) : x(_x), y(_y), z(_z), w(_w), r(_r), g(_g), b(_b), a(_a) {}

    // xy + rgb + depth, for vertex converted to the screen
    vertex(double _x, double _y, double _r, double _g, double _b, double _depth) : x(_x), y(_y), r(_r), g(_g), b(_b), depth(_depth) {}

    // xy + rgb, for pixel on the screen
    vertex(double _x, double _y, double _r, double _g, double _b) : x(_x), y(_y), r(_r), g(_g), b(_b) {}

    void sRGB_to_linear_RGB()
    {
        r = sRGB_to_linear(r);
        g = sRGB_to_linear(g);
        b = sRGB_to_linear(b);
    }

    void linear_RGB_to_SRGB()
    {
        r = linear_to_SRGB(r);
        g = linear_to_SRGB(g);
        b = linear_to_SRGB(b);
    }

private:
    int sRGB_to_linear(int channel)
    {
        double linear = channel / 255.0;
        if (linear <= 0.04045)
        {
            linear = linear / 12.92;
        }
        else
        {
            linear = pow((linear + 0.055) / 1.055, 2.4);
        }
        return (int)(linear * 255.0 + 0.5);
    }

    int linear_to_SRGB(int channel)
    {
        double sRGB = channel / 255.0;
        if (sRGB <= 0.0031308)
        {
            sRGB = sRGB * 12.92;
        }
        else
        {
            sRGB = 1.055 * pow(sRGB, 1.0 / 2.4) - 0.055;
        }
        return (int)(sRGB * 255.0 + 0.5);
    }
};

// read input files
vector<vector<string>> read_file(const string &fileName)
{
    ifstream input;
    input.open(fileName);
    vector<vector<string>> contents;

    if (input.is_open())
    {
        string aline;
        while (getline(input, aline))
        {
            vector<string> content;
            string word = "";
            size_t i = 0;

            // remove the space in the front
            while (i < aline.size() && (aline[i] == ' ' || aline[i] == '\t'))
            {
                i++;
            }

            for (; i < aline.size(); i++)
            {
                if (aline[i] == ' ' || aline[i] == '\t')
                {
                    while (i < aline.size() && (aline[i] == ' ' || aline[i] == '\t'))
                    {
                        i++;
                    }
                    i--;
                    content.push_back(word);
                    word = "";
                }
                else
                {
                    word = word + aline[i];
                }
            }

            if (word != "")
            {
                content.push_back(word);
            }

            if (content.size() >= 1)
            {
                // skip the empty line
                contents.push_back(content);
            }
        }
        input.close();
    }
    else
    {
        cerr << "Can't open file: " << fileName << "." << endl;
        exit(-1);
    }
    return contents;
}

// compare 2 vertices in y axis
bool compare_vertex(vertex a, vertex b)
{
    return (a.y < b.y);
}

// swap 2 vertices
void swap_vertex(vertex &a, vertex &b)
{
    vertex tmp(b);
    b.x = a.x;
    b.y = a.y;
    b.r = a.r;
    b.g = a.g;
    b.b = a.b;
    b.depth = a.depth;
    b.a = a.a;

    a.x = tmp.x;
    a.y = tmp.y;
    a.r = tmp.r;
    a.g = tmp.g;
    a.b = tmp.b;
    a.depth = tmp.depth;
    a.a = tmp.a;
}

vertex operator+(const vertex &a, const vertex &b)
{
    vertex result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    result.w = a.w + b.w;
    result.r = a.r + b.r;
    result.g = a.g + b.g;
    result.b = a.b + b.b;
    result.depth = a.depth + b.depth;
    result.a = a.a + b.a;

    return result;
}

vertex operator-(const vertex &a, const vertex &b)
{
    vertex result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    result.w = a.w - b.w;
    result.r = a.r - b.r;
    result.g = a.g - b.g;
    result.b = a.b - b.b;
    result.depth = a.depth - b.depth;
    result.a = a.a - b.a;

    return result;
}

vertex operator*(const vertex &a, const double c)
{
    vertex result;
    result.x = a.x * c;
    result.y = a.y * c;
    result.z = a.z * c;
    result.w = a.w * c;
    result.r = a.r * c;
    result.g = a.g * c;
    result.b = a.b * c;
    result.depth = a.depth * c;
    result.a = a.a * c;

    return result;
}

vertex operator/(const vertex &a, const double c)
{
    vertex result;
    result.x = a.x / c;
    result.y = a.y / c;
    result.z = a.z / c;
    result.w = a.w / c;
    result.r = a.r / c;
    result.g = a.g / c;
    result.b = a.b / c;
    result.depth = a.depth / c;
    result.a = a.a / c;

    return result;
}

inline int sign(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

void draw_line(vertex p0, vertex p1, vector<vertex> &pixels, int width, int height)
{

    double xDiff = p1.x - p0.x;
    double yDiff = p1.y - p0.y;
    double steps = std::max(fabs(xDiff), fabs(yDiff));

    double xStep = xDiff / steps;
    double yStep = yDiff / steps;

    double x = p0.x + 0.5 * sign(xStep);
    double y = p0.y + 0.5 * sign(yStep);

    int lastX = round(p0.x), lastY = round(p0.y);
    int currX, currY;

    double rStep = (p1.r - p0.r) / steps;
    double gStep = (p1.g - p0.g) / steps;
    double bStep = (p1.b - p0.b) / steps;
    double r = p0.r, g = p0.g, b = p0.b;

    for (int i = 0; i <= steps; i++)
    {
        currX = round(x);
        currY = round(y);

        if (currX >= 0 && currX < width && currY >= 0 && currY < height)
        {
            if (currX != lastX || currY != lastY)
            {
                // double dist = sqrt((currX - x) * (currX - x) + (currY - y) * (currY - y));
                // vertex new_pix(lastX, lastY, r + dist * rStep, g + dist * gStep, b + dist * bStep);
                vertex new_pix(lastX, lastY, r, g, b);
                pixels.push_back(new_pix);
            }

            lastX = currX;
            lastY = currY;
        }

        lastX = currX;
        lastY = currY;

        x += xStep;
        y += yStep;

        r += rStep;
        g += gStep;
        b += bStep;
    }
}

void draw_point(const vertex p0, double pointsize, vector<vertex> &pixels, double **depth_map)
{
    double x0 = p0.x;
    double y0 = p0.y;

    double half_width = pointsize / 2;
    int xmin = ceil(x0 - half_width);
    int xmax = ceil(x0 + half_width);
    int ymin = ceil(y0 - half_width);
    int ymax = ceil(y0 + half_width);

    for (int x = xmin; x < xmax; ++x)
    {
        for (int y = ymin; y < ymax; ++y)
        {
            if (x >= 0 && x < width && y >= 0 && y < height)
            {
                vertex p(x, y, p0.r, p0.g, p0.b);

                // check depth
                if (depth_map[x][y] > p0.depth && p0.depth >= -1)
                {
                    depth_map[x][y] = p0.depth;
                    pixels.push_back(p);
                }
            }
        }
    }
}

void blend(vertex &new_pixel, vertex *orignal)
{
    double color_s[4] = {new_pixel.r, new_pixel.g, new_pixel.b, new_pixel.a};
    double color_d[4] = {orignal->r, orignal->g, orignal->b, orignal->a};

    color_s[3] /= 255;
    color_d[3] /= 255;
    // a' = a_s + a_d * (1 - a_s)
    double a_prime = color_s[3] + color_d[3] * (1 - color_s[3]);
    new_pixel.r = (color_s[3] / a_prime) * color_s[0] + ((1 - color_s[3]) * color_d[3] / a_prime) * color_d[0];
    new_pixel.g = (color_s[3] / a_prime) * color_s[1] + ((1 - color_s[3]) * color_d[3] / a_prime) * color_d[1];
    new_pixel.b = (color_s[3] / a_prime) * color_s[2] + ((1 - color_s[3]) * color_d[3] / a_prime) * color_d[2];
    new_pixel.a = a_prime * 255;
}

bool clip_line(const vertex p0, const vertex p1, vertex &clipped_p0, vertex &clipped_p1)
{
    // input p0, p1
    bool is_p0 = false;
    bool is_p1 = false;
    // if ((p0.x <= p0.w) && (p0.x >= -p0.w) && (p0.y <= p0.w) && (p0.y >= -p0.w) && (p0.z <= p0.w) && (p0.z >= -p0.w))
    if (p0.w > p0.x)
    {
        // p0 is in the view
        is_p0 = true;
        clipped_p0 = p0;
        clipped_p0.x = (clipped_p0.x / clipped_p0.w + 1) * (width / 2.0);
        clipped_p0.y = (clipped_p0.y / clipped_p0.w + 1) * (height / 2.0);
    }

    // if ((p1.x <= p1.w) && (p1.x >= -p1.w) && (p1.y <= p1.w) && (p1.y >= -p1.w) && (p1.z <= p1.w) && (p1.z >= -p1.w))
    if (p1.w > p1.x)
    {
        // p1 is in the view
        is_p1 = true;
        clipped_p1 = p1;
        clipped_p1.x = (clipped_p1.x / clipped_p1.w + 1) * (width / 2.0);
        clipped_p1.y = (clipped_p1.y / clipped_p1.w + 1) * (height / 2.0);
    }

    if (is_p0 && is_p1)
    {
        return true;
    }

    // if (!is_p0 && !is_p1)
    // {
    //     // the line is out of the screen
    //     return false;
    // }

    double dist_0;
    double dist_1;

    if (!is_p0)
    // while (!(p0.x <= p0.w) && (p0.x >= -p0.w) && (p0.y <= p0.w) && (p0.y >= -p0.w) && (p0.z <= p0.w) && (p0.z >= -p0.w))
    {
        // p0 is off-screen
        dist_0 = p0.x - p0.w;
        dist_1 = p1.x - p1.w;

        clipped_p0 = p1 * (dist_0 / (dist_0 - dist_1)) - p0 * (dist_1 / (dist_0 - dist_1));
        clipped_p1 = p1;
    }

    if (!is_p1)
    {
        dist_0 = p0.x - p0.w;
        dist_1 = p1.x - p1.w;

        clipped_p0 = p0;
        clipped_p1 = p0 * (dist_1 / (dist_1 - dist_0)) - p1 * (dist_0 / (dist_1 - dist_0));

        // p0 = clipped_p0;
        // p1 = clipped_p1;
    }

    clipped_p0.x = (clipped_p0.x / clipped_p0.w + 1) * (width / 2.0);
    clipped_p0.y = (clipped_p0.y / clipped_p0.w + 1) * (height / 2.0);
    clipped_p1.x = (clipped_p1.x / clipped_p1.w + 1) * (width / 2.0);
    clipped_p1.y = (clipped_p1.y / clipped_p1.w + 1) * (height / 2.0);

    if (clipped_p0.y > clipped_p1.y)
        swap(clipped_p0, clipped_p1);

    return true;
}

bool if_counter_clock(vertex v1, vertex v2, vertex v3)
{
    // Check if the triangle is back-facing
    double edge1_x = v2.x - v1.x;
    double edge1_y = v2.y - v1.y;
    double edge2_x = v3.x - v1.x;
    double edge2_y = v3.y - v1.y;
    double cross_z = edge1_x * edge2_y - edge1_y * edge2_x;

    if (cross_z <= 0)
        return false;

    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cerr << "Usage: ./program file.txt" << endl;
        return -1;
    }

    auto contents = read_file(argv[1]);

    // image information
    string filename;

    bool is_depth = false;
    bool is_sRGB = false;
    bool is_rgba = false;
    bool is_frustum = false;
    bool is_cull = false;
    // bool is_hyp = false;

    // rgba buffer
    map<pair<double, double>, vertex *> rgb_buffer;

    // store all the vertices
    vector<vertex> vertices;

    // all the triangles need to draw, each triangle containts 3 vertices
    vector<vector<vertex>> triangles;

    // draw line
    vector<vector<vertex>> lines;

    // draw point
    vector<vertex> point_center;
    vector<double> point_size;

    // store all the pixels for drawing
    vector<vertex> pixels;

    // color to fit in
    double current_color[4] = {255.0, 255.0, 255.0, 255.0};

    // scan input
    for (size_t j = 0; j < contents.size(); j++)
    {
        // for each line
        if (contents[j].size() >= 4 && contents[j][0] == "png")
        {
            width = stoi(contents[j][1]);
            height = stoi(contents[j][2]);
            filename = contents[j][3];
        }
        else if (contents[j].size() >= 5 && contents[j][0] == "xyzw")
        {
            vertex new_vertex(stod(contents[j][1]), stod(contents[j][2]), stod(contents[j][3]), stod(contents[j][4]), current_color[0], current_color[1], current_color[2], current_color[3]);
            vertices.push_back(new_vertex);
        }
        else if (contents[j].size() >= 4 && contents[j][0] == "tri")
        {
            vector<vertex> new_tri;
            for (size_t i = 1; i < contents[j].size(); i++)
            {
                // get the vertex
                int vertex_num = stoi(contents[j][i]);
                if (vertex_num > 0)
                {
                    vertex_num--;
                }
                else
                {
                    vertex_num = vertices.size() + vertex_num;
                }

                double x = vertices[vertex_num].x;
                double y = vertices[vertex_num].y;
                double w = vertices[vertex_num].w;
                double z = vertices[vertex_num].z;

                vertex new_vertex((x / w + 1.0) * (width / 2.0), (y / w + 1.0) * (height / 2.0), vertices[vertex_num].r, vertices[vertex_num].g, vertices[vertex_num].b, z / w);

                if (is_sRGB)
                {
                    // sRGB to RGB
                    new_vertex.sRGB_to_linear_RGB();
                }
                if (is_rgba)
                {
                    new_vertex.a = vertices[vertex_num].a;
                }
                if (is_frustum)
                {
                    new_tri.push_back(vertices[vertex_num]);
                }
                else
                {
                    new_tri.push_back(new_vertex);
                }
            }

            if (is_frustum)
            {
                // move it to tri part
                sort(new_tri.begin(), new_tri.end(), compare_vertex);

                vector<vertex> frustum_triangle;
                vector<vertex> frustum_triangle_1;
                vector<vertex> frustum_triangle_2;
                vector<vertex> frustum_triangle_3;
                vector<vertex> frustum_triangle_4;

                vertex p0 = new_tri[0];
                vertex p1 = new_tri[1];
                vertex p2 = new_tri[2];

                if (p0.w < 0 && p1.w < 0 && p2.w < 0)
                {
                    // triangle is out of screen
                    // do nothing
                }
                else
                {
                    vertex clipped_p0;
                    vertex clipped_p1;
                    vertex clipped_p2;
                    vertex clipped_p3;
                    vertex clipped_p4;
                    vertex clipped_p5;

                    bool clip_p0_p1 = clip_line(p0, p1, clipped_p0, clipped_p1);
                    bool clip_p0_p2 = clip_line(p0, p2, clipped_p2, clipped_p3);
                    bool clip_p1_p2 = clip_line(p1, p2, clipped_p4, clipped_p5);

                    if (clip_p0_p1 && !clip_p0_p2 && !clip_p1_p2)
                    {
                        frustum_triangle.empty();
                        frustum_triangle.push_back(p0);
                        frustum_triangle.push_back(p1);
                        frustum_triangle.push_back(clipped_p5);
                        triangles.push_back(frustum_triangle);

                        frustum_triangle.empty();
                        frustum_triangle.push_back(p0);
                        frustum_triangle.push_back(clipped_p3);
                        frustum_triangle.push_back(clipped_p5);
                        triangles.push_back(frustum_triangle);
                    }
                    if (clip_p0_p2 && !clip_p1_p2 && !clip_p0_p1)
                    {
                        frustum_triangle.empty();
                        frustum_triangle.push_back(p0);
                        frustum_triangle.push_back(clipped_p4);
                        frustum_triangle.push_back(clipped_p1);
                        triangles.push_back(frustum_triangle);

                        frustum_triangle.empty();
                        frustum_triangle.push_back(p0);
                        frustum_triangle.push_back(p2);
                        frustum_triangle.push_back(clipped_p4);
                        triangles.push_back(frustum_triangle);
                    }
                    if (clip_p1_p2 && !clip_p0_p2 && !clip_p0_p1)
                    {
                        frustum_triangle.empty();
                        frustum_triangle.push_back(p1);
                        frustum_triangle.push_back(clipped_p0);
                        frustum_triangle.push_back(clipped_p2);
                        triangles.push_back(frustum_triangle);

                        frustum_triangle.empty();
                        frustum_triangle.push_back(p1);
                        frustum_triangle.push_back(p2);
                        frustum_triangle.push_back(clipped_p2);
                        triangles.push_back(frustum_triangle);
                    }
                    // if (!clip_p0_p2)
                    // {
                    //     frustum_triangle.empty();
                    //     frustum_triangle.push_back(clipped_p0);
                    //     frustum_triangle.push_back(clipped_p5);
                    //     frustum_triangle.push_back(p1);
                    //     triangles.push_back(frustum_triangle);
                    // }
                    // if (!clip_p1_p2)
                    // {
                    //     frustum_triangle.empty();
                    //     frustum_triangle.push_back(clipped_p1);
                    //     frustum_triangle.push_back(clipped_p3);
                    //     frustum_triangle.push_back(p0);
                    //     triangles.push_back(frustum_triangle);
                    // }
                    else
                    {
                        frustum_triangle_1.push_back(clipped_p0);
                        frustum_triangle_1.push_back(clipped_p1);
                        frustum_triangle_1.push_back(clipped_p4);

                        frustum_triangle_2.push_back(clipped_p0);
                        frustum_triangle_2.push_back(clipped_p2);
                        frustum_triangle_2.push_back(clipped_p4);

                        frustum_triangle_3.push_back(clipped_p2);
                        frustum_triangle_3.push_back(clipped_p5);
                        frustum_triangle_3.push_back(clipped_p4);

                        frustum_triangle_4.push_back(clipped_p2);
                        frustum_triangle_4.push_back(clipped_p5);
                        frustum_triangle_4.push_back(clipped_p3);

                        triangles.push_back(frustum_triangle_1);
                        triangles.push_back(frustum_triangle_2);
                        triangles.push_back(frustum_triangle_3);
                        triangles.push_back(frustum_triangle_4);
                    }
                }
            }
            else
            {
                if (is_cull)
                {
                    if (!if_counter_clock(new_tri[0], new_tri[1], new_tri[2]))
                    {
                        triangles.push_back(new_tri);
                    }
                }
                else
                {
                    triangles.push_back(new_tri);
                }
            }
        }
        else if (contents[j].size() >= 4 && contents[j][0] == "rgb")
        {
            current_color[0] = 1.0 * stoi(contents[j][1]);
            current_color[1] = 1.0 * stoi(contents[j][2]);
            current_color[2] = 1.0 * stoi(contents[j][3]);
        }
        else if (contents[j].size() >= 4 && contents[j][0] == "rgba")
        {
            is_rgba = true;
            current_color[0] = 1.0 * stoi(contents[j][1]);
            current_color[1] = 1.0 * stoi(contents[j][2]);
            current_color[2] = 1.0 * stoi(contents[j][3]);
            current_color[3] = 255.0 * stod(contents[j][4]);
        }
        else if (contents[j].size() >= 1 && contents[j][0] == "depth")
        {
            is_depth = true;
        }
        else if (contents[j].size() >= 1 && contents[j][0] == "sRGB")
        {
            is_sRGB = true;
        }
        else if (contents[j].size() >= 1 && contents[j][0] == "frustum")
        {
            is_frustum = true;
        }
        else if (contents[j].size() >= 1 && contents[j][0] == "cull")
        {
            is_cull = true;
        }
        else if (contents[j].size() >= 3 && contents[j][0] == "line")
        {
            vector<vertex> new_line;
            for (size_t i = 1; i < contents[j].size(); i++)
            {
                // get the vertex
                int vertex_num = stoi(contents[j][i]);
                if (vertex_num > 0)
                {
                    vertex_num--;
                }
                else
                {
                    vertex_num = vertices.size() + vertex_num;
                }

                double x = vertices[vertex_num].x;
                double y = vertices[vertex_num].y;
                double w = vertices[vertex_num].w;
                double z = vertices[vertex_num].z;

                vertex point((x / w + 1.0) * (width / 2.0), (y / w + 1.0) * (height / 2.0), vertices[vertex_num].r, vertices[vertex_num].g, vertices[vertex_num].b, z / w);

                new_line.push_back(point);
            }
            lines.push_back(new_line);
        }
        else if (contents[j].size() >= 3 && contents[j][0] == "point")
        {
            point_size.push_back(stod(contents[j][1]));
            // get the vertex
            int vertex_num = stoi(contents[j][2]);
            if (vertex_num > 0)
            {
                vertex_num--;
            }
            else
            {
                vertex_num = vertices.size() + vertex_num;
            }

            double x = vertices[vertex_num].x;
            double y = vertices[vertex_num].y;
            double w = vertices[vertex_num].w;
            double z = vertices[vertex_num].z;

            vertex point((x / w + 1.0) * (width / 2.0), (y / w + 1.0) * (height / 2.0), vertices[vertex_num].r, vertices[vertex_num].g, vertices[vertex_num].b, z / w);

            point_center.push_back(point);
        }
    }

    // create a 2D arrary for depth mapping
    double **depth_map = new double *[height];
    for (int i = 0; i < height; i++)
    {
        depth_map[i] = new double[width];
    }

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            depth_map[i][j] = 1;
        }
    }

    // for each triangles, use DDA algorithm
    for (size_t i = 0; i < triangles.size(); i++)
    {
        // sort vertex
        vector<vertex> triangle(triangles[i]);
        sort(triangle.begin(), triangle.end(), compare_vertex);

        // DDA algorithm
        vertex p0 = triangle[0];
        vertex p1 = triangle[1];
        vertex p2 = triangle[2];

        // vector a, b, c
        vertex a = p1 - p0;
        vertex b = p2 - p0;
        vertex c = p2 - p1;

        // if there's no triangle
        if (a.x == b.x && a.y == b.y)
        {
            continue;
        }
        else
        {
            // for p0-p1 -> p0-p2
            if (a.x > b.x)
            {
                swap_vertex(a, b);
            }

            // calculate step
            vertex step_a = a / a.y;
            vertex step_b = b / b.y;

            double init_y = ceil(p0.y);
            double delta_y = init_y - p0.y;

            // init point for the vert first line's start and end
            vertex init_a = p0 + step_a * delta_y;
            vertex init_b = p0 + step_b * delta_y;

            int y = int(init_y);
            // scan line
            for (int k = 0; int(init_y) + k < p1.y; k++)
            {
                y = int(init_y) + k;

                // start and the point of the line
                vertex start = init_a + step_a * k;
                vertex end = init_b + step_b * k;

                int start_x = max(int(ceil(start.x)), 0);
                int end_x = min(int(ceil(end.x)), width);

                // if (end_x == int(end_x_double))
                // {
                //     end_x++;
                // }

                vertex delta = (end - start) / (end.x - start.x);

                for (int x = start_x; x < end_x; x++)
                {
                    vertex new_pixel = start + delta * (x - start.x);

                    // if (is_hyp)
                    // {
                    //     double w = new_pixel.w;
                    //     new_pixel.r /= w;
                    //     new_pixel.g /= w;
                    //     new_pixel.b /= w;
                    //     new_pixel.w = 1 / w;
                    // }
                    if (is_sRGB && !is_rgba)
                    {
                        // RGB to sRGB
                        new_pixel.linear_RGB_to_SRGB();
                        // pixels.push_back(new_pixel);
                    }

                    if (is_depth)
                    {
                        // check depth
                        if (depth_map[x][y] > new_pixel.depth && new_pixel.depth >= -1)
                        {
                            depth_map[x][y] = new_pixel.depth;
                            pixels.push_back(new_pixel);
                        }
                    }
                    else if (is_rgba)
                    {
                        // caculate new rgb
                        if (rgb_buffer.count(make_pair(new_pixel.x, new_pixel.y)) == 1)
                        {
                            vertex *orignal = rgb_buffer[make_pair(new_pixel.x, new_pixel.y)];
                            blend(new_pixel, orignal);
                        }
                        rgb_buffer[make_pair(new_pixel.x, new_pixel.y)] = &new_pixel;
                        int flag = 0;
                        for (vector<vertex>::iterator it = pixels.begin(); it != pixels.end(); ++it) // const auto pix : pixels)
                        {
                            if (it->x == new_pixel.x && it->y == new_pixel.y)
                            {
                                flag = 1;
                                it->r = new_pixel.r;
                                it->g = new_pixel.g;
                                it->b = new_pixel.b;
                                it->a = new_pixel.a;
                            }
                        }
                        if (flag == 0)
                        {
                            pixels.push_back(new_pixel);
                        }
                    }
                    else
                    {
                        pixels.push_back(new_pixel);
                    }
                }
            }

            // for p1-p2 -> p0->p2
            // b may be swapped
            b = p2 - p0;

            // point on line p0->p2 where y = p1
            vertex pm = (p2 - p0) / (p2.y - p0.y) * (p1.y - p0.y) + p0;

            // make sure p1,c go first
            if (pm.x < p1.x)
            {
                swap_vertex(c, b);
                swap_vertex(pm, p1);
            }

            vertex step_c = c / c.y;
            step_b = b / b.y;

            init_y = ceil(p1.y);
            delta_y = init_y - p1.y;

            vertex init_c = p1 + step_c * delta_y;
            init_b = pm + step_b * delta_y;

            y = int(init_y);
            // scan line
            for (int k = 0; int(init_y) + k < p2.y; k++)
            {
                y = int(init_y) + k;

                // start and the point of the line
                vertex start = init_c + step_c * k;
                vertex end = init_b + step_b * k;

                int start_x = max(int(ceil(start.x)), 0);
                int end_x = min(int(ceil(end.x)), width);

                // if (end_x == int(end_x_double))
                // {
                //     end_x++;
                // }

                vertex delta = (end - start) / (end.x - start.x);

                for (int x = start_x; x < end_x; x++)
                {
                    vertex new_pixel = start + delta * (x - start.x);

                    if (is_sRGB && !is_rgba)
                    {
                        new_pixel.linear_RGB_to_SRGB();
                    }
                    if (is_depth)
                    {
                        // check depth
                        if (depth_map[x][y] > new_pixel.depth && new_pixel.depth >= -1)
                        {
                            depth_map[x][y] = new_pixel.depth;
                            pixels.push_back(new_pixel);
                        }
                    }
                    else if (is_rgba)
                    {
                        // new rgba
                        if (rgb_buffer.count(make_pair(new_pixel.x, new_pixel.y)) == 1)
                        {
                            vertex *orignal = rgb_buffer[make_pair(new_pixel.x, new_pixel.y)];
                            blend(new_pixel, orignal);
                        }
                        // push to map
                        rgb_buffer[make_pair(new_pixel.x, new_pixel.y)] = &new_pixel;

                        // cout << new_pixel.r << ", " << new_pixel.g << ", " << new_pixel.b << ", " << new_pixel.a << endl;

                        int flag = 0;
                        for (size_t i = 0; i < pixels.size(); i++)
                        {
                            if (pixels[i].x == new_pixel.x && pixels[i].y == new_pixel.y)
                            {
                                flag = 1;
                                pixels[i] = new_pixel;
                            }
                        }
                        if (flag == 0)
                        {
                            pixels.push_back(new_pixel);
                        }
                    }
                    else
                    {
                        pixels.push_back(new_pixel);
                    }
                }
            }
        }
    }

    if (is_rgba)
    {
        for (size_t i = 0; i < pixels.size(); i++)
        {
            pixels[i].linear_RGB_to_SRGB();
        }
    }

    // line
    for (size_t i = 0; i < lines.size(); i++)
    {
        // for lines
        draw_line(lines[i][0], lines[i][1], pixels, width, height);
    }

    // point
    for (size_t i = 0; i < point_size.size(); i++)
    {
        // void draw_point(const vertex &p0, double pointsize, vector<vertex> &pixels, double **depth_map)
        draw_point(point_center[i], point_size[i], pixels, depth_map);
    }
    // generate images
    unsigned char *image = (unsigned char *)calloc(width * height * 4, sizeof(unsigned char));
    for (size_t i = 0; i < pixels.size(); i++)
    {
        int x = pixels[i].x;
        int y = pixels[i].y;
        int red = pixels[i].r;
        int green = pixels[i].g;
        int blue = pixels[i].b;
        int alpha = pixels[i].a;

        image[((y * width) + x) * 4 + 0] = red;
        image[((y * width) + x) * 4 + 1] = green;
        image[((y * width) + x) * 4 + 2] = blue;
        image[((y * width) + x) * 4 + 3] = alpha;
    }

    lodepng_encode32_file(filename.c_str(), image, width, height);
    free(image);

    // Deallocate memory
    for (int i = 0; i < height; i++)
    {
        delete[] depth_map[i];
    }
    delete[] depth_map;

    return 0;
}