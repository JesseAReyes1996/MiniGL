/*
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <stack>

using namespace std;

/*
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

struct vertex
{
    vec3 color;
    vec4 pos;
};

struct triangle
{
    vertex a, b, c;
};

MGLpoly_mode draw_mode;

vec3 current_color;
vector<triangle> LoT;
vector<vertex> LoV;

MGLmatrix_mode mat_mode;
mat4 projectionMatrix;
mat4 modelviewMatrix;

stack<mat4> pmStack;
stack<mat4> mvmStack;

vector<vector<MGLfloat>> zBuffer;
/*
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description)
{
    printf("%s\n", description);
    exit(1);
}

/*
 * Calculate the area of a triangle.
 */
MGLfloat areaOfATriangle(triangle current)
{
    return current.a.pos[0] * ((current.b.pos[1] - current.c.pos[1])) + (current.a.pos[1] * (current.c.pos[0] - current.b.pos[0])) + ((current.b.pos[0] * current.c.pos[1]) - (current.b.pos[1] * current.c.pos[0]));
}

/*
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{

    int red;
    int green;
    int blue;

    MGLfloat dx = static_cast<float>(width) / 2.0; // pixel width
    MGLfloat dy = static_cast<float>(height) / 2.0;// pixel height

    MGLfloat x;
    MGLfloat y;

    MGLfloat i;
    MGLfloat j;

    MGLfloat areaABC;

    MGLfloat alpha;
    MGLfloat beta;
    MGLfloat gamma;

    vertex pixelCoord;
    vector<vertex> pixelCoords;

    triangle ABC;
    triangle PPP;
    vertex P;

    vector<MGLfloat> zVal;
    MGLfloat zDepth;
    //Resize the 2D buffer
    zBuffer.resize(width);
    for(uint i = 0; i < width; ++i)
    {
        zBuffer[i].resize(height,INFINITY);
    }

    for(uint iter = 0; iter < LoT.size(); ++iter)
    {
        x = LoT.at(iter).a.pos[0] / LoT.at(iter).a.pos[3];
        y = LoT.at(iter).a.pos[1] / LoT.at(iter).a.pos[3];
        zVal.push_back(LoT.at(iter).a.pos[2] / LoT.at(iter).a.pos[3]);

        //Pixel coordinates A
        i = (x + 1.0) * dx - 0.5;
        j = (y + 1.0) * dy - 0.5;

        pixelCoord.pos = vec4(i, j, 0, 1);
        pixelCoords.push_back(pixelCoord);

        x = LoT.at(iter).b.pos[0] / LoT.at(iter).b.pos[3];
        y = LoT.at(iter).b.pos[1] / LoT.at(iter).b.pos[3];
        zVal.push_back(LoT.at(iter).b.pos[2] / LoT.at(iter).b.pos[3]);

        //Pixel coordinates B
        i = ((x + 1.0) * dx) - 0.5;
        j = ((y + 1.0) * dy) - 0.5;

        pixelCoord.pos = vec4(i, j, 0, 1);
        pixelCoords.push_back(pixelCoord);

        x = LoT.at(iter).c.pos[0] / LoT.at(iter).c.pos[3];
        y = LoT.at(iter).c.pos[1] / LoT.at(iter).c.pos[3];
        zVal.push_back(LoT.at(iter).c.pos[2] / LoT.at(iter).c.pos[3]);

        //Pixel coordinates C
        i = ((x + 1.0) * dx) - 0.5;
        j = ((y + 1.0) * dy) - 0.5;

        pixelCoord.pos = vec4(i, j, 0, 1);
        pixelCoords.push_back(pixelCoord);

        ABC.a = pixelCoords.at(pixelCoords.size() - 3);
        ABC.b = pixelCoords.at(pixelCoords.size() - 2);
        ABC.c = pixelCoords.at(pixelCoords.size() - 1);

        areaABC = areaOfATriangle(ABC);

        //Define a bounding box per triangle
        int minX = max(min(min(pixelCoords.at(pixelCoords.size() - 3).pos[0], pixelCoords.at(pixelCoords.size() - 2).pos[0]), pixelCoords.at(pixelCoords.size() - 1).pos[0]), static_cast<float>(0));
        int minY = max(min(min(pixelCoords.at(pixelCoords.size() - 3).pos[1], pixelCoords.at(pixelCoords.size() - 2).pos[1]), pixelCoords.at(pixelCoords.size() - 1).pos[1]), static_cast<float>(0));
        int maxX = min(max(max(pixelCoords.at(pixelCoords.size() - 3).pos[0], pixelCoords.at(pixelCoords.size() - 2).pos[0]), pixelCoords.at(pixelCoords.size() - 1).pos[0]), static_cast<float>(width));
        int maxY = min(max(max(pixelCoords.at(pixelCoords.size() - 3).pos[1], pixelCoords.at(pixelCoords.size() - 2).pos[1]), pixelCoords.at(pixelCoords.size() - 1).pos[1]), static_cast<float>(height));

        //Calculate the barycentric coordinates of the triangle and if
        //a point lies in/on the triangle, fill it in.
        //for(uint h = 0; h < height; ++h)
        //{
        //    for(uint w = 0; w < width; ++w)
        //    {

        for(int h = minY; h < maxY; ++h)
        {
            for(int w = minX; w < maxX; ++w)
            {
                P.pos = vec4(w, h, 0, 1);

                PPP.a = P;
                PPP.b = ABC.b;
                PPP.c = ABC.c;
                alpha = areaOfATriangle(PPP) / areaABC;

                PPP.a = ABC.a;
                PPP.b = P;
                PPP.c = ABC.c;
                beta  = areaOfATriangle(PPP) / areaABC;

                PPP.a = ABC.a;
                PPP.b = ABC.b;
                PPP.c = P;
                gamma = areaOfATriangle(PPP) / areaABC;

                if((alpha >= 0 && alpha <= 1) && (beta >= 0 && beta <= 1) && (gamma >= 0 && gamma <= 1))
                {
                    MGLfloat tAlpha = (alpha/LoT.at(iter).a.pos[3])/((alpha/LoT.at(iter).a.pos[3])+(beta/LoT.at(iter).b.pos[3])+(gamma/LoT.at(iter).c.pos[3]));
                    MGLfloat tBeta  = (beta/LoT.at(iter).b.pos[3])/((alpha/LoT.at(iter).a.pos[3])+(beta/LoT.at(iter).b.pos[3])+(gamma/LoT.at(iter).c.pos[3]));
                    MGLfloat tGamma = (gamma/LoT.at(iter).c.pos[3])/((alpha/LoT.at(iter).a.pos[3])+(beta/LoT.at(iter).b.pos[3])+(gamma/LoT.at(iter).c.pos[3]));

                    zDepth = tAlpha*zVal.at(zVal.size() - 3) + tBeta*zVal.at(zVal.size() - 2) + tGamma*zVal.at(zVal.size() - 1);

                    red   = tAlpha*round(LoT.at(iter).a.color[0]) + tBeta*round(LoT.at(iter).b.color[0]) + tGamma*round(LoT.at(iter).c.color[0]);
                    green = tAlpha*round(LoT.at(iter).a.color[1]) + tBeta*round(LoT.at(iter).b.color[1]) + tGamma*round(LoT.at(iter).c.color[1]);
                    blue  = tAlpha*round(LoT.at(iter).a.color[2]) + tBeta*round(LoT.at(iter).b.color[2]) + tGamma*round(LoT.at(iter).c.color[2]);

                    if(zDepth < zBuffer.at(w).at(h) && zDepth >= -1 && zDepth <= 1)
                    {
                        zBuffer.at(w).at(h) = zDepth;
                        data[w + h * width] = Make_Pixel(red, green, blue);
                    }
                }
            }
        }
    }
    zBuffer.clear();
}

/*
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode (MGL_TRIANGLES or MGL_QUADS).
 */
void mglBegin(MGLpoly_mode mode)
{
    draw_mode = mode;
}

/*
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    triangle temp;

    if(draw_mode == MGL_TRIANGLES)
    {
        for(uint i = 0; i+2 < LoV.size(); i += 3)
        {
            temp.a = LoV.at(i);
            temp.b = LoV.at(i+1);
            temp.c = LoV.at(i+2);
            LoT.push_back(temp);
        }
    }

    if(draw_mode == MGL_QUADS)
    {
        for(uint i = 0; i+3 < LoV.size(); i += 4)
        {
            temp.a = LoV.at(i);   //A
            temp.b = LoV.at(i+1); //B
            temp.c = LoV.at(i+2); //C
            LoT.push_back(temp);

            temp.a = LoV.at(i);   //A
            temp.b = LoV.at(i+2); //C
            temp.c = LoV.at(i+3); //D
            LoT.push_back(temp);
        }
    }
    LoV.clear();
}

/*
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x, y, 0);
}

/*
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vertex v;
    v.color = current_color;
    v.pos = vec4(x, y, z, 1);

    //Multiply the projection matrix by the modelview matrix by the vertex
    //before appending it to the list of vertices to apply any transformations
    //stored in the projection matrix and the modelview matrix.
    v.pos = projectionMatrix * modelviewMatrix * v.pos;

    LoV.push_back(v);
}

/*
 * Set the current matrix mode (MGL_MODELVIEW or MGL_PROJECTION).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    mat_mode = mode;
}

/*
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    if(mat_mode == MGL_PROJECTION)
    {
        pmStack.push(projectionMatrix);
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        mvmStack.push(modelviewMatrix);
    }
}

/*
 * Assign the top of the stack for the current matrix mode to the current matrix
 * and then pop the top of the stack for the current matrix mode.
 */
void mglPopMatrix()
{
    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = pmStack.top();
        pmStack.pop();
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = mvmStack.top();
        mvmStack.pop();
    }
}

/*
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix.make_zero();
        projectionMatrix.values[0] = 1;
        projectionMatrix.values[5] = 1;
        projectionMatrix.values[10] = 1;
        projectionMatrix.values[15] = 1;
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix.make_zero();
        modelviewMatrix.values[0] = 1;
        modelviewMatrix.values[5] = 1;
        modelviewMatrix.values[10] = 1;
        modelviewMatrix.values[15] = 1;
    }
}

/*
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    if(mat_mode == MGL_PROJECTION)
    {
        for(uint i = 0; i < 16; ++i)
        {
            projectionMatrix.values[i] = matrix[i];
        }
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        for(uint i = 0; i < 16; ++i)
        {
            modelviewMatrix.values[i] = matrix[i];
        }
    }
}

/*
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
    mat4 multMat;
    multMat.make_zero();

    for(uint i = 0; i < 16; ++i)
    {
        multMat.values[i] = matrix[i];
    }

    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = projectionMatrix * multMat;

    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = modelviewMatrix * multMat;
    }
}

/*
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mat4 translateMatrix;
    translateMatrix.make_zero();

    translateMatrix.values[0]  = 1;
    translateMatrix.values[5]  = 1;
    translateMatrix.values[10] = 1;
    translateMatrix.values[12] = x;
    translateMatrix.values[13] = y;
    translateMatrix.values[14] = z;
    translateMatrix.values[15] = 1;

    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = projectionMatrix * translateMatrix;
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = modelviewMatrix * translateMatrix;
    }
}

/*
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    mat4 rotationMatrix;
    rotationMatrix.make_zero();

    //Convert degrees to radians
    MGLfloat radians = (M_PI/180) * angle;

    //Normalize the vector
    MGLfloat magnitude = sqrt((x*x) + (y*y) + (z*z));
    x = x / magnitude;
    y = y / magnitude;
    z = z / magnitude;

    rotationMatrix.values[0]  = (x*x)*(1-cos(radians)) + cos(radians);
    rotationMatrix.values[1]  = (y*x)*(1-cos(radians)) + z*sin(radians);
    rotationMatrix.values[2]  = (z*x)*(1-cos(radians)) - y*sin(radians);
    rotationMatrix.values[4]  = (x*y)*(1-cos(radians)) - z*sin(radians);
    rotationMatrix.values[5]  = (y*y)*(1-cos(radians)) + cos(radians);
    rotationMatrix.values[6]  = (z*y)*(1-cos(radians)) + x*sin(radians);
    rotationMatrix.values[8]  = (x*z)*(1-cos(radians)) + y*sin(radians);
    rotationMatrix.values[9]  = (y*z)*(1-cos(radians)) - x*sin(radians);
    rotationMatrix.values[10] = (z*z)*(1-cos(radians)) + cos(radians);
    rotationMatrix.values[15] = 1;

    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = projectionMatrix * rotationMatrix;
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = modelviewMatrix * rotationMatrix;
    }
}

/*
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mat4 scaleMatrix;
    scaleMatrix.make_zero();

    scaleMatrix.values[0]  = x;
    scaleMatrix.values[5]  = y;
    scaleMatrix.values[10] = z;
    scaleMatrix.values[15] = 1;

    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = projectionMatrix * scaleMatrix;
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = modelviewMatrix * scaleMatrix;
    }
}

/*
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
    mat4 perspectiveMatrix;
    perspectiveMatrix.make_zero();

    perspectiveMatrix.values[0]  = ((2.0 * near)/(right - left));
    perspectiveMatrix.values[5]  = ((2.0 * near)/(top - bottom));
    perspectiveMatrix.values[8]  = (right + left)/(right - left);
    perspectiveMatrix.values[9]  = (top + bottom)/(top - bottom);
    perspectiveMatrix.values[10] = -((far + near)/(far - near));
    perspectiveMatrix.values[11] = -1.0;
    perspectiveMatrix.values[14] = -((2 * far * near)/(far - near));

    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = projectionMatrix * perspectiveMatrix;
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = modelviewMatrix * perspectiveMatrix;
    }
}

/*
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    mat4 orthoMatrix;
    orthoMatrix.make_zero();

    orthoMatrix.values[0]  = ((2.0)/(right - left));
    orthoMatrix.values[5]  = ((2.0)/(top - bottom));
    orthoMatrix.values[10] = ((-2.0)/(far - near));
    orthoMatrix.values[12] = -((right + left)/(right - left));
    orthoMatrix.values[13] = -((top + bottom)/(top - bottom));
    orthoMatrix.values[14] = -((far + near)/(far - near));
    orthoMatrix.values[15] = 1.0;

    if(mat_mode == MGL_PROJECTION)
    {
        projectionMatrix = projectionMatrix * orthoMatrix;
    }

    if(mat_mode == MGL_MODELVIEW)
    {
        modelviewMatrix = modelviewMatrix * orthoMatrix;
    }
}

/*
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
    current_color = vec3(round(red*255), round(green*255), round(blue*255));
}
