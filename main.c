#include <SDL2/SDL.h>
#include <stdio.h>
#include <math.h>

#define Cw 600
#define Ch 600

struct light
{
    double *position;
    double *direction;
    double intensity;
};

double intersectRaySphere(double O[3], double V[3], double C[3], double r);
double dotProduct(double v1[3], double v2[3]);
void vectorPointSubtract(double *resultVector, double coordinate1[3], double coordinate2[3]);
void vectorPointAdd(double *resultVector, double coordinate1[3], double coordinate2[3]);
void traceRay(int color[4], double O[3], double V[3], double C[3], double radius, struct light lights[], int numLights, double specular, double reflective);
void canvasToViewport(double *V, int x, int y);
double lightReflected(double N[3], double lightDir[3]);
int clamp(int color[3]);

int main(int argc, char *argv[])
{
    // declare and initialize
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;

    // initialize the video library
    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        printf("video initialize error: %s\n", SDL_GetError());
        return 1;
    }

    // set the window and renderer
    if (SDL_CreateWindowAndRenderer(Cw, Ch, 0, &window, &renderer) == -1)
    {
        printf("render initialize error: %s\n", SDL_GetError());
        return 1;
    }

    // draw points
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    SDL_RenderClear(renderer);

    double origin[3] = {0, -0.3, 0.35};
    double center[3] = {0, 0, 3};
    double radius = 1;
    double specular = 30;
    double reflective = 0.8;

    struct light lights[3];
    lights[0].intensity = 0.2;
    lights[0].position = NULL;
    lights[0].direction = NULL;

    lights[1].intensity = 0.6;
    lights[1].position = (double *)malloc(3 * sizeof(double));
    lights[1].position[0] = 2, lights[1].position[1] = 1, lights[1].position[2] = 0;
    lights[1].direction = NULL;

    lights[2].intensity = 0.2;
    lights[2].position = NULL;
    lights[2].direction = (double *)malloc(3 * sizeof(double));
    lights[2].direction[0] = 1, lights[2].direction[1] = 4, lights[2].direction[2] = 4;

    for (int x = -Cw / 2; x <= Cw / 2; x++)
    {
        for (int y = -Ch / 2; y <= Ch; y++)
        {
            double viewport[3];

            canvasToViewport(viewport, x, y);

            int color[4];
            traceRay(color, origin, viewport, center, radius, lights, sizeof(lights) / sizeof(lights[0]), specular, reflective);

            SDL_SetRenderDrawColor(renderer, color[0], color[1], color[2], color[4]);
            SDL_RenderDrawPoint(renderer, Cw / 2 + x, Ch / 2 - y);
        }
    }

    // present the renderer on the window
    SDL_RenderPresent(renderer);
    SDL_Delay(30000);
    //}

    // end process
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    free(lights[1].position);
    free(lights[2].direction);

    return 0;
}

/// @brief solves for t in the equation P = O + t(V - O), the parameter/scalar that determines where on the line the surface point lies
/// @param O the position of the camera (considered as the origin of the rendered scene)
/// @param V The position of the viewport, which is the 3d canvas that represents the 2d pixel canvas
/// @param C The position of the center of the sphere
/// @param r The radius of the sphere
/// @return t value
double intersectRaySphere(double O[3], double V[3], double C[3], double r)
{

    // P = O + t(V - O)
    // D = (V - O)
    double D[3];
    vectorPointSubtract(D, V, O);

    // r = P - C  =>  r^2 = ||P - C||^2
    // r^2 = ||O + tD - C||^2
    // CO = O - C
    double CO[3];
    vectorPointSubtract(CO, O, C);

    // (||D||^2)t^2 + 2||CO*D||t + ||CO||^2 - r^2
    double a = dotProduct(D, D);
    double b = 2 * dotProduct(CO, D);
    double c = dotProduct(CO, CO) - r * r;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0)
    {
        return -1;
    }
    else
    {
        double t1 = (-b + sqrt(discriminant)) / (2 * a);
        double t2 = (-b - sqrt(discriminant)) / (2 * a);

        if (t2 > 0)
        {
            return t2;
        }
        else if (t1 > 0)
        {
            return t1;
        }
        else
        {
            return -1;
        }
    }
}

/// @brief calculates dot product of two vectors
/// @param v1 vector 1
/// @param v2 vector 2
/// @return dot product
double dotProduct(double v1[3], double v2[3])
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/// @brief Subtracts two vectors or two points, since the calculation is identical
/// @param resultVector The resultant vector from the subtraction, represented as an array
/// @param coordinate1 vector 1 or point 1
/// @param coordinate2 vector 2 or point 2
void vectorPointSubtract(double *resultVector, double coordinate1[3], double coordinate2[3])
{
    resultVector[0] = coordinate1[0] - coordinate2[0];
    resultVector[1] = coordinate1[1] - coordinate2[1];
    resultVector[2] = coordinate1[2] - coordinate2[2];
}

/// @brief Adds two vectors or two points, since the calculation is identical
/// @param resultVector The resultant vector from the addition, represented as an array
/// @param coordinate1 vector 1 or point 1
/// @param coordinate2 vector 2 or point 2
void vectorPointAdd(double *resultVector, double coordinate1[3], double coordinate2[3])
{
    resultVector[0] = coordinate1[0] + coordinate2[0];
    resultVector[1] = coordinate1[1] + coordinate2[1];
    resultVector[2] = coordinate1[2] + coordinate2[2];
}

/// @brief Based on the given t value (distance from the camera to a surface point), returns the corresponding color value
/// @param color array that contains the rgb values of the color at the point
/// @param t parameter that determines how far along the line (that extends from the camera and beyond the canvas) the point is
void traceRay(int color[4], double O[3], double V[3], double C[3], double r, struct light lights[], int numLights, double specular, double reflective)
{
    double t = intersectRaySphere(O, V, C, r);

    double D[3];
    vectorPointSubtract(D, V, O);
    double P[3];
    P[0] = O[0] + t * D[0];
    P[1] = O[1] + t * D[1];
    P[2] = O[2] + t * D[2];

    double floorHeight = -1;
    double tFloor = 0;

    if (D[1] != 0)
    {
        tFloor = (floorHeight - O[1]) / D[1];
    }

    if (t != -1 && P[1] > floorHeight)
    {
        // normal vector
        double N[3];
        vectorPointSubtract(N, P, C);
        double lengthN = sqrt(dotProduct(N, N));
        N[0] = N[0] / lengthN;
        N[1] = N[1] / lengthN;
        N[2] = N[2] / lengthN;

        double totalIntensity = 0;

        for (int i = 0; i < numLights; i++)
        {
            int s = specular;

            // dot product of normal and direction vector
            double L[3] = {0, 0, 0};

            if (lights[i].position != NULL)
            {
                vectorPointSubtract(L, lights[i].position, P);
                totalIntensity += lights[i].intensity * lightReflected(N, L);
            }
            else if (lights[i].direction != NULL)
            {
                vectorPointAdd(L, L, lights[i].direction);
                totalIntensity += lights[i].intensity * lightReflected(N, L);
            }
            else
            {
                totalIntensity += lights[i].intensity;
            }

            if (s != -1)
            {
                double R[3];
                double n_dot_l = dotProduct(N, L);
                R[0] = 2 * n_dot_l * N[0] - L[0];
                R[1] = 2 * n_dot_l * N[1] - L[1];
                R[2] = 2 * n_dot_l * N[2] - L[2];

                double view[3];
                view[0] = -1 * D[0];
                view[1] = -1 * D[1];
                view[2] = -1 * D[2];

                double lr = lightReflected(R, view);
                double power = pow(lr, specular);
                totalIntensity += lights[i].intensity * power;
            }
        }

        if (reflective <= 0)
        {
            color[0] = 160 * totalIntensity;
            color[1] = 160 * totalIntensity;
            color[2] = 160 * totalIntensity;
            color[3] = 255;
        }
        else
        {
            double oppositeD[3] = {0, 0, 0};
            vectorPointSubtract(oppositeD, oppositeD, D);

            double n_dot_minusD = dotProduct(N, oppositeD);

            double reflectedRay[3];
            reflectedRay[0] = 2 * n_dot_minusD * N[0] - oppositeD[0];
            reflectedRay[1] = 2 * n_dot_minusD * N[1] - oppositeD[1];
            reflectedRay[2] = 2 * n_dot_minusD * N[2] - oppositeD[2];

            int reflectedColor[4];
            vectorPointAdd(reflectedRay, reflectedRay, P);

            traceRay(reflectedColor, P, reflectedRay, C, 0, lights, numLights, 0, 0);

            color[0] = 160 * totalIntensity * (1 - reflective) + (reflectedColor[0] * reflective);
            color[1] = 160 * totalIntensity * (1 - reflective) + (reflectedColor[1] * reflective);
            color[2] = 160 * totalIntensity * (1 - reflective) + (reflectedColor[2] * reflective);
            color[3] = 255;
        }

        clamp(color);
    }
    else if (tFloor > 0)
    {
        double shadowIntensity = 0;

        double floorPoint[3];
        floorPoint[0] = O[0] + tFloor * D[0];
        floorPoint[1] = O[1] + tFloor * D[1];
        floorPoint[2] = O[2] + tFloor * D[2];

        for (int i = 0; i < numLights; i++)
        {
            double tFloorSphere = 0;

            if (lights[i].position == NULL && lights[i].direction == NULL)
            {
                shadowIntensity += lights[i].intensity;
            }
            else
            {
                if (lights[i].position != NULL)
                {
                    tFloorSphere = intersectRaySphere(floorPoint, lights[i].position, C, r);
                }
                else
                {
                    double floorToLight[3] = {0, 0, 0};
                    vectorPointSubtract(floorToLight, floorToLight, lights[i].direction);
                    vectorPointAdd(floorToLight, floorToLight, floorPoint);
                    tFloorSphere = intersectRaySphere(floorPoint, floorToLight, C, r);
                }

                if (tFloorSphere == -1 || tFloorSphere > 1)
                {
                    shadowIntensity += lights[i].intensity;
                }
            }
        }

        if ((int)floorPoint[0] % 2 == 0)
        {
            color[0] = 50 * shadowIntensity;
            color[1] = 50 * shadowIntensity;
            color[2] = 50 * shadowIntensity;
            color[3] = 255;
        }
        else
        {
            color[0] = 255 * shadowIntensity;
            color[1] = 255 * shadowIntensity;
            color[2] = 255 * shadowIntensity;
            color[3] = 255;
        }
    }
    else
    {
        color[0] = 174;
        color[1] = 198;
        color[2] = 207;
        color[3] = 255;
    }
}

void canvasToViewport(double *V, int x, int y)
{
    double Vw = 1;
    double Vh = 1;
    double z = 1;

    V[0] = x * Vw / Cw;
    V[1] = y * Vh / Ch;
    V[2] = z;
}

double lightReflected(double N[3], double lightDir[3])
{
    double dot = dotProduct(lightDir, N);
    if (dot > 0)
    {
        return dot / (sqrt(dotProduct(lightDir, lightDir)) * sqrt(dotProduct(N, N)));
    }
    else
    {
        return 0;
    }
}

int clamp(int color[3])
{
    for (int i = 0; i <= 3; i++)
    {
        if (color[i] < 0)
        {
            color[i] = 0;
        }
        if (color[i] > 255)
        {
            color[i] = 255;
        }
    }
}