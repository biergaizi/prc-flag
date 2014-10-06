/*
* prc_flag.c - draw the national flag of PRC
*
* $ cc prc_flag.c -I/usr/include/cairo -lcairo -lm -std=c99
*
* Copyright (C) 2014 Tom Li. All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#if __STDC_VERSION__ > 199900L
#define __IS_C99_COMPILER
#endif
#ifndef __IS_C99_COMPILER
#error "Please use a C99 compiler (or add '-std=c99') to compile"
#endif

#ifndef _GNU_SOURCE /* M_PI is a (GNU) extension */
#define _GNU_SOURCE /* It was disabled by default in C99 */
#endif

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <cairo.h>


typedef struct Point {
    double x;
    double y;
} Point;

#define POINT(a, b) \
(Point) \
{ \
  .x = (a), \
  .y = (b) \
}

#define WIDTH 960
#define HEIGHT 640

#define BIG_STAR POINT(5, 5)

#define BIG_STAR_R 3
#define SMALL_STAR_R 1

#define SMALL_STAR_1 POINT(10, 2)
#define SMALL_STAR_2 POINT(12, 4)
#define SMALL_STAR_3 POINT(12, 7)
#define SMALL_STAR_4 POINT(10, 9)

#define BASE (HEIGHT / 2 / 10)


void star(cairo_t *cr, Point *center, double r, Point *circle);
void rotate(cairo_t *cr, Point *point, double radian);
Point circle_line_intersection(Point *circle, double r, Point *p1, Point *p2);
double radian_two_points(Point *p1, Point *p2, double r);
void grid(cairo_t *cr, Point *topleft, Point *botright, unsigned int row, unsigned int col);
void line(cairo_t *cr, Point *p1, Point *p2);
void flag_star(cairo_t *cr, Point *p, double r);
void flag_line(cairo_t *cr, Point *center);


void star(cairo_t *cr, Point *center, double r, Point *circle)
{
    Point topleft   = POINT(center->x - r * sin(2 * M_PI / 5),
                            center->y - r * cos(2 * M_PI / 5));

    Point topright  = POINT(center->x + r * sin(2 * M_PI / 5),
                            center->y - r * cos(2 * M_PI / 5));

    Point botleft   = POINT(center->x - r * sin(M_PI / 5),
                            center->y + r * cos(M_PI / 5));

    Point topcenter = POINT(center->x,
                            center->y - r);

    Point botright  = POINT(center->x + r * sin(M_PI / 5),
                            center->y + r * cos(M_PI / 5));

    cairo_save(cr);
    if (circle) {
        Point point = circle_line_intersection(center, r, circle, center);
        double rad = radian_two_points(&topleft, &point, r);
        rotate(cr, center, rad);
    }
    cairo_new_path(cr);
    cairo_move_to(cr, topleft.x, topleft.y);
    cairo_line_to(cr, topright.x, topright.y);
    cairo_line_to(cr, botleft.x, botleft.y);
    cairo_line_to(cr, topcenter.x, topcenter.y);
    cairo_line_to(cr, botright.x, botright.y);
    cairo_close_path(cr);
    cairo_set_source_rgb(cr, 255, 255, 0);
    cairo_fill(cr);
    cairo_restore(cr);
}

void rotate(cairo_t *cr, Point *point, double radian)
{
    cairo_translate(cr, point->x, point->y);
    cairo_rotate(cr, radian);
    cairo_translate(cr, -point->x, -point->y);
}

Point circle_line_intersection(Point *circle, double r, Point *p1, Point *p2)
{
    /* y = mx + c, solve m and c from p1 and p2 */
    double lm = (p2->y - p1->y) / (p2->x - p1->x);
    double lc = p1->y - p1->x * lm;
    assert(lc == p2->y - p2->x * lm);

    /* solve (x - p) ^ 2 + (y - q) ^ 2 = r ^ 2 and y = mx + c */
    double a = pow(lm, 2) + 1;
    double b = 2 * ((lm * lc) - (lm * circle->y) - circle->x);
    double c = pow(circle->y, 2) - pow(r, 2) + pow(circle->x, 2) - (2 * lc * circle->y) + pow(lc, 2);

    double delta = pow(b, 2) - 4 * a * c;
    assert(delta >= 0);
    delta = sqrt(delta);

    double x1 = (-b - delta) / (2 * a);
    double y1 = lm * x1 + lc;

    return POINT(x1, y1);
}

double radian_two_points(Point *p1, Point *p2, double r)
{
    double distance = sqrt(pow(p2->x - p1->x, 2) + pow(p2->y - p1->y, 2));
    double angle = 2 * asin(distance / 2 / r);
    if (angle > M_PI && p2->y < p1->y) {
        angle = -angle;
    }
    else if (angle < M_PI && p2->y > p1->y) {
        angle = -angle;
    }
    return angle;
}

void grid(cairo_t *cr, Point *topleft, Point *botright, unsigned int row, unsigned int col)
{
    cairo_save(cr);
    cairo_set_source_rgb(cr, 0, 0, 0);  /* black */
    cairo_set_line_width(cr, 2);

    float max_y = (botright->y - topleft->y);
    float max_x  = (botright->x - topleft->x);

    for (int i = 0; i < max_y; i += max_y / row) {
        cairo_move_to(cr, i, topleft->y);
        cairo_line_to(cr, i, botright->x);
    }

    for (int i = 0; i < max_x; i += max_x / col) {
        cairo_move_to(cr, topleft->x, i);
        cairo_line_to(cr, botright->y, i);
    }
    cairo_stroke(cr);
    cairo_restore(cr);
}

void line(cairo_t *cr, Point *p1, Point *p2)
{
    cairo_save(cr);
    cairo_set_source_rgb(cr, 0, 0, 0);  /* black */
    cairo_move_to(cr, p1->x, p1->y);
    cairo_line_to(cr, p2->x, p2->y);
    cairo_stroke(cr);
    cairo_restore(cr);
}

void flag_star(cairo_t *cr, Point *p, double r)
{
    Point *circle = NULL;
    if (p->x != BIG_STAR.y && p->y != BIG_STAR.y) {
        circle = &POINT(BIG_STAR.x * BASE,
                        BIG_STAR.y * BASE);
        circle->x = BIG_STAR.x * BASE;
        circle->y = BIG_STAR.y * BASE;
    }
    p->x *= BASE;
    p->y *= BASE;
    r *= BASE;

    star(cr, p, r, circle);
}

void flag_line(cairo_t *cr, Point *center)
{
    line(cr,
         &POINT(BIG_STAR.x * BASE,
                BIG_STAR.y * BASE),
         &POINT(center->x * BASE,
                center->y * BASE));
}

void show_usage(char *name)
{
    printf("Usage: %s [OPTION]\n", name);
    printf("   -a\tDraw all auxiliary lines\n");
    printf("   -o [FILE]\tOutput file (will auto-append \".png\" if not ends with it)\n");
    printf("   -h\tThis help message\n");
}

int str_ends_with(const char *str, const char *suffix)
{
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix > lenstr) {
        return 0;
    }
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

void *xmalloc(size_t size)
{
    static char error_msg[] = "malloc() failed\n";
    void *ptr = malloc(size);
    if (!ptr) {
        fputs(error_msg, stderr);
        abort();
    }
    return ptr;
}

int main(int argc, char **argv)
{
    static struct option longopts[] = {
        {"auxiliary-lines", no_argument,       0, 'a'},
        {"output",          required_argument, 0, 'o'},
        {"help",            no_argument,       0, 'h'},
    };

    bool auxiliary_line = false;
    char *output_arg = NULL;
    char *output_real = NULL;

    int opt;
    while ((opt = getopt_long(argc, argv, "hao:", longopts, NULL)) != -1) {
        switch (opt) {
            case 'a':
                auxiliary_line = true;
                break;
            case 'o':
                output_arg = optarg;
                break;
            case '?':
            case 'h':
                show_usage(argv[0]);
                exit(0);
        }
    }
    if (!output_arg) {
        show_usage(argv[0]);
        exit(0);
    }
    else {
        if (!str_ends_with(output_arg, ".png")) {
            output_real = xmalloc(strlen(output_arg) + strlen(".png"));
            sprintf(output_real, "%s%s", output_arg, ".png");
            printf("renamed '%s' to '%s'\n", output_arg, output_real);
        }
        else {
            output_real = xmalloc(strlen(output_arg));
            strcpy(output_real, output_arg);
        }
    }

    cairo_surface_t *surface;
    cairo_t *cr;

    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, WIDTH, HEIGHT);
    cr = cairo_create(surface);

    /* fill background to red */
    cairo_set_source_rgb(cr, 255, 0, 0);  /* red */
    cairo_rectangle(cr, 0, 0, WIDTH, HEIGHT);
    cairo_fill(cr);

    flag_star(cr, &BIG_STAR, BIG_STAR_R);
    flag_star(cr, &SMALL_STAR_1, SMALL_STAR_R);
    flag_star(cr, &SMALL_STAR_2, SMALL_STAR_R);
    flag_star(cr, &SMALL_STAR_3, SMALL_STAR_R);
    flag_star(cr, &SMALL_STAR_4, SMALL_STAR_R);

    if (auxiliary_line) {
        flag_line(cr, &SMALL_STAR_1);
        flag_line(cr, &SMALL_STAR_2);
        flag_line(cr, &SMALL_STAR_3);
        flag_line(cr, &SMALL_STAR_4);

        grid(cr, &POINT(0, 0), &POINT(HEIGHT, WIDTH),
             2, 2);
        grid(cr, &POINT(0, 0), &POINT(HEIGHT / 2, WIDTH / 2),
             15, 10);
    }

    assert(output_real != NULL);
    cairo_status_t status = cairo_surface_write_to_png(surface, output_real);
    if (status != CAIRO_STATUS_SUCCESS) {
        fprintf(stderr, "cairo_surface_write_to_png() failed: %s\n", cairo_status_to_string(status));
        exit(1);
    }

    free(output_real);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

    return 0;
}
