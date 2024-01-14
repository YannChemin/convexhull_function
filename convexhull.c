#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define SPECTRALBANDSMAX 200

// Structure to represent a point in 2D space
typedef struct {
    double x;
    double y;
} Point;


// Function to calculate the orientation of three points (p, q, r)
double orientation(Point p, Point q, Point r) {
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
}

// Function to find the convex hull of a set of points
void convexHull(Point *points, int n, Point *hull, int *hullSize) {
    // Find the leftmost point
    int l = 0;

    // Find the point with maximum x-coordinate
    int maxXIndex = 0;
    for (int i = 1; i < n; ++i) {
        if (points[i].x > points[maxXIndex].x) {
            maxXIndex = i;
        }
    }

    // Storage for the convex hull
    *hullSize = 0;

    // Start from the leftmost point and move clockwise until the point with maximum x-coordinate
    int p = l;
    int q;
    while (p != maxXIndex) {
        // Add the current point to the result
        hull[(*hullSize)++] = points[p];

        // Search for a point 'q' such that orientation(p, q, x) is clockwise for all points 'x'
        q = (p + 1) % n;
        for (int i = 0; i < n; ++i) {
            if (orientation(points[p], points[i], points[q]) > 0) {
                q = i;
            }
        }

        // Now q is the most clockwise with respect to p
        // Set p as q for the next iteration
        p = q;
    }

    // Add the point with maximum x-coordinate to the result
    hull[(*hullSize)++] = points[maxXIndex];
}

// Comparison function for qsort
int comparePoints(const void *a, const void *b) {
    return ((Point *)a)->x > ((Point *)b)->x ? 1 : -1;
}

// Function to perform linear interpolation
double interpolate(double x, Point p1, Point p2) {
    double denominator = p2.x - p1.x;
    if (denominator == 0.0) {
        // Avoid division by zero, return a default value (e.g., p1.y)
        if (p1.y > 0) return p1.y;
        else if (p2.y > 0) return p2.y;
        else return 0;
    } else { 
        return p1.y + (x - p1.x) * (p2.y - p1.y) / denominator;
    }
}


// Function to interpolate and fill y-coordinates
void interpolateAndFill(Point *originalData, int originalSize, Point *convexHull, int hullSize) {
    // Sort convex hull points based on x-coordinates
    qsort(convexHull, hullSize, sizeof(Point), comparePoints);

    // Interpolate linearly to fill y-coordinates of original data
    int j = 0;
    for (int i = 0; i < originalSize; ++i) {
        double x = originalData[i].x;

        // Find the convex hull points within the x-range [x1, x2]
        while (j < hullSize - 1 && convexHull[j + 1].x <= x) {
            j++;
        }

        // Interpolate and update the y-coordinate
        double outval = interpolate(x, convexHull[j], convexHull[j + 1]);
        if (outval > 0){
            originalData[i].y = outval;
        }
    }
}

void continuumremoval(Point *points, int numPoints, Point *cr){
    // Make a backup for later
    Point orig[SPECTRALBANDSMAX];
    for (int i = 0; i < numPoints; ++i) {
        orig[i].x = points[i].x;
        orig[i].y = points[i].y;
    }
    // CONVEX HULL COMPUTATION
    // Create arrays to store the convex hull points
    Point hull[numPoints];
    int hullSize;
    // Call the convexHull function
    convexHull(points, numPoints, hull, &hullSize);
    // INTERPOLATION OF CONVEX HULL OUTPUT TO INPUT X-COORD LENGTH
    // Call the interpolateAndFill function
    interpolateAndFill(points, numPoints, hull, hullSize);
    // Removal of continuum = original / convex hull
    for (int i = 0; i < numPoints; ++i) {
        cr[i].x = orig[i].x;
        cr[i].y = orig[i].y / points[i].y;
    }
}

int main() {
    // Read data into an array of points
    Point cr[SPECTRALBANDSMAX];
    Point points[SPECTRALBANDSMAX];

    int numPoints = 19;
    //wavelengths,values
    points[0].x=2.1017;
    points[0].y=40.356634632034634;
    points[1].x=2.1197;
    points[1].y=40.573279653679656;
    points[2].x=2.1374;
    points[2].y=40.529950649350646;
    points[3].x=2.1548;
    points[3].y=40.685935064935066;
    points[4].x=2.172;
    points[4].y=40.183318614718615;
    points[5].x=2.1893;
    points[5].y=39.31673852813853;
    points[6].x=2.2068;
    points[6].y=37.982205194805196;
    points[7].x=2.2254;
    points[7].y=37.33227012987013;
    points[8].x=2.2417;
    points[8].y=36.812322077922076;
    points[9].x=2.2592;
    points[9].y=36.51768484848485;
    points[10].x=2.2756;
    points[10].y=36.11905800865801;
    points[11].x=2.2928;
    points[11].y=35.63377316017316;
    points[12].x=2.3097;
    points[12].y=34.74119567099567;
    points[13].x=2.326;
    points[13].y=33.70129956709957;
    points[14].x=2.3421;
    points[14].y=33.138022510822516;
    points[15].x=2.3588;
    points[15].y=32.68740086580087;
    points[16].x=2.3745;
    points[16].y=32.91271168831169;
    points[17].x=2.3903;
    points[17].y=33.11202510822511;
    points[18].x=2.4064;
    points[18].y=33.24201212121212;

    int i;
    #pragma omp parallel for default(none) private(i) shared(points, numPoints, cr)
    for (i = 0; i < 1000000; i++){
        continuumremoval(points, numPoints, cr);
    }

    // Write the continuum removed points to a file in gnuplot .dat format
    FILE *outFile = fopen("cr_output.dat", "w");
    if (outFile == NULL) {
        perror("Error creating output file");
        return 1;
    }
    for (i = 0; i < numPoints; ++i) {
        fprintf(outFile, "%lf\t%lf\n", cr[i].x, cr[i].y);
    }
    // Close the output file
    fclose(outFile);

    return 0;
}

