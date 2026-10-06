#ifndef LIGHTHELIXFIT_H
#define LIGHTHELIXFIT_H

#include <array>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include "Hit.h"

namespace LightHelixFit
{
struct Point3D
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// This mirrors the small type aliases in TrackFitUtils.h L19-L23.
using position_t = std::pair<double, double>;
using position_vector_t = std::vector<position_t>;
using circle_fit_output_t = std::tuple<double, double, double>; // [R, x0, y0]
using line_fit_output_t = std::tuple<double, double>;            // [slope, intercept]
using fitpars_t = std::array<double, 5>;                         // [R, X0, Y0, zslope, Z0]

struct FitResult
{
    bool success = false;

    // Same parameter ordering as TrackFitUtils::fitClusters:
    // [R, X0, Y0, zslope, Z0]. See TrackFitUtils.h L128-L130 and
    // TrackFitUtils.cc L646-L650.
    fitpars_t fitpars = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};

    position_vector_t rz_points;
    position_t circle_pca_xy = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    Point3D helix_pca = {};
    double dca2d = std::numeric_limits<double>::quiet_NaN();
    double dca3d = std::numeric_limits<double>::quiet_NaN();
};

static inline Point3D MakePoint3D(Hit *hit)
{
    if (!hit)
        return {};
    return {hit->posX(), hit->posY(), hit->posZ()};
}

static inline double square(const double x) { return x * x; }

static inline double RadiusXY(const Point3D &p)
{
    return std::hypot(p.x, p.y);
}

static inline double Distance3D(const Point3D &a, const Point3D &b)
{
    return std::sqrt(square(a.x - b.x) + square(a.y - b.y) + square(a.z - b.z));
}

static inline position_vector_t MakeRZPoints(const Point3D &vertex, const Point3D &hit1, const Point3D &hit2)
{
    // This corresponds to the R-Z input preparation used by TrackFitUtils::line_fit
    // in TrackFitUtils.cc L272-L280, except here the three points are fixed to
    // {vertex, hit1, hit2}.
    position_vector_t rz_points;
    rz_points.reserve(3);
    rz_points.emplace_back(RadiusXY(vertex), vertex.z);
    rz_points.emplace_back(RadiusXY(hit1), hit1.z);
    rz_points.emplace_back(RadiusXY(hit2), hit2.z);
    return rz_points;
}

static inline circle_fit_output_t circle_fit_by_three_points(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    // This is the lightweight replacement for the circle fit stage declared in
    // TrackFitUtils.h L22-L40 and used in TrackFitUtils.cc L622.
    // TrackFitUtils uses Taubin for many-point fits; here we solve the exact
    // three-point circle because the intended input is only {vertex, hit1, hit2}.
    const double x1 = p0.x;
    const double y1 = p0.y;
    const double x2 = p1.x;
    const double y2 = p1.y;
    const double x3 = p2.x;
    const double y3 = p2.y;

    const double det = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    if (std::abs(det) < 1e-9)
    {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    const double r1sq = square(x1) + square(y1);
    const double r2sq = square(x2) + square(y2);
    const double r3sq = square(x3) + square(y3);

    const double x0 = (r1sq * (y2 - y3) + r2sq * (y3 - y1) + r3sq * (y1 - y2)) / det;
    const double y0 = (r1sq * (x3 - x2) + r2sq * (x1 - x3) + r3sq * (x2 - x1)) / det;
    const double R = std::hypot(x1 - x0, y1 - y0);

    return {R, x0, y0};
}

static inline line_fit_output_t line_fit(const position_vector_t &positions)
{
    // This follows the same Deming-style line fit implementation used in
    // TrackFitUtils.cc L238-L268.
    double xmean = 0.0;
    double ymean = 0.0;
    for (const auto &[x, y] : positions)
    {
        xmean += x;
        ymean += y;
    }

    const double n = positions.size();
    if (n <= 0.0)
    {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    xmean /= n;
    ymean /= n;

    double ssd_x = 0.0;
    double ssd_y = 0.0;
    double ssd_xy = 0.0;
    for (const auto &[x, y] : positions)
    {
        ssd_x += square(x - xmean);
        ssd_y += square(y - ymean);
        ssd_xy += (x - xmean) * (y - ymean);
    }

    if (std::abs(ssd_xy) < 1e-12)
    {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    const double slope = (ssd_y - ssd_x + std::sqrt(square(ssd_y - ssd_x) + 4.0 * square(ssd_xy))) / (2.0 * ssd_xy);
    const double intercept = ymean - slope * xmean;
    return {slope, intercept};
}

static inline position_t get_circle_point_pca(const double radius, const double x0, const double y0, const Point3D &global)
{
    // This corresponds directly to TrackFitUtils.h L126 and TrackFitUtils.cc L596-L605.
    const position_t origin(x0, y0);
    const position_t point(global.x, global.y);

    const double dx = point.first - origin.first;
    const double dy = point.second - origin.second;
    const double norm = std::hypot(dx, dy);
    if (norm <= 0.0 || !std::isfinite(norm) || !std::isfinite(radius))
    {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }

    return {origin.first + radius * dx / norm,
            origin.second + radius * dy / norm};
}

static inline double get_helix_pathlength(const fitpars_t &fitpars, const Point3D &start_point, const Point3D &end_point)
{
    // This mirrors the path-length helper in TrackFitUtils.h L152-L153 and
    // TrackFitUtils.cc L905-L937.
    const double R = fitpars[0];
    const double x0 = fitpars[1];
    const double y0 = fitpars[2];
    const double zslope = fitpars[3];

    const double half_turn_z_pathlength = zslope * 2.0 * R;
    const double z_dist = end_point.z - start_point.z;

    int n_turns = 0;
    if (std::abs(half_turn_z_pathlength) > 1e-12)
        n_turns = static_cast<int>(std::floor(std::abs(z_dist / half_turn_z_pathlength)));

    const double phic_start = std::atan2(start_point.y - y0, start_point.x - x0);
    const double phic_end = std::atan2(end_point.y - y0, end_point.x - x0);
    const double phic_dist = std::abs(phic_end - phic_start) + n_turns * 2.0 * M_PI;
    const double xy_dist = R * phic_dist;

    return std::sqrt(xy_dist * xy_dist + z_dist * z_dist);
}

static inline Point3D get_helix_pca(const fitpars_t &fitpars, const Point3D &global)
{
    // This follows the same PCA construction used in TrackFitUtils.h L124-L126
    // and TrackFitUtils.cc L559-L593.
    const double radius = fitpars[0];
    const double x0 = fitpars[1];
    const double y0 = fitpars[2];
    const double zslope = fitpars[3];
    const double z0 = fitpars[4];

    const position_t pca_circle = get_circle_point_pca(radius, x0, y0, global);
    if (!std::isfinite(pca_circle.first) || !std::isfinite(pca_circle.second))
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    const double pca_circle_radius = std::hypot(pca_circle.first, pca_circle.second);
    const double pca_z = pca_circle_radius * zslope + z0;
    return {pca_circle.first, pca_circle.second, pca_z};
}

static inline FitResult fitClusters(const Point3D &vertex, const Point3D &hit1, const Point3D &hit2)
{
    // This is the lightweight analogue of TrackFitUtils::fitClusters declared in
    // TrackFitUtils.h L128-L130 and implemented in TrackFitUtils.cc L611-L652.
    // The structural difference is intentional: TrackFitUtils works with an arbitrary
    // cluster vector plus cluskeys, while this helper is specialized to exactly
    // {vertex, INTT hit 1, INTT hit 2} and therefore skips the cluskey/use_intt logic.
    FitResult out;

    const circle_fit_output_t circle_fit_pars = circle_fit_by_three_points(vertex, hit1, hit2);
    const double R = std::get<0>(circle_fit_pars);
    const double X0 = std::get<1>(circle_fit_pars);
    const double Y0 = std::get<2>(circle_fit_pars);
    if (!std::isfinite(R) || !std::isfinite(X0) || !std::isfinite(Y0) || R <= 0.0)
        return out;

    out.rz_points = MakeRZPoints(vertex, hit1, hit2);
    const line_fit_output_t line_fit_pars = line_fit(out.rz_points);
    const double zslope = std::get<0>(line_fit_pars);
    const double Z0 = std::get<1>(line_fit_pars);
    if (!std::isfinite(zslope) || !std::isfinite(Z0))
        return out;

    out.fitpars = {R, X0, Y0, zslope, Z0};
    out.circle_pca_xy = get_circle_point_pca(R, X0, Y0, vertex);
    out.helix_pca = get_helix_pca(out.fitpars, vertex);
    if (!std::isfinite(out.circle_pca_xy.first) || !std::isfinite(out.circle_pca_xy.second) || !std::isfinite(out.helix_pca.z))
        return out;

    out.dca2d = std::hypot(out.circle_pca_xy.first - vertex.x, out.circle_pca_xy.second - vertex.y);
    out.dca3d = Distance3D(out.helix_pca, vertex);
    out.success = std::isfinite(out.dca2d) && std::isfinite(out.dca3d);
    return out;
}

static inline FitResult fitClusters(const Point3D &vertex, Hit *hit1, Hit *hit2)
{
    if (!hit1 || !hit2)
        return {};
    return fitClusters(vertex, MakePoint3D(hit1), MakePoint3D(hit2));
}

// Convenience wrappers kept for analysis readability. These are not named like the
// TrackFitUtils entry points, but they forward to the fitClusters-style interface above.
static inline FitResult FitFromThreePoints(const Point3D &vertex, const Point3D &hit1, const Point3D &hit2)
{
    return fitClusters(vertex, hit1, hit2);
}

static inline FitResult FitFromVertexAndDoublet(const Point3D &vertex, Hit *hit1, Hit *hit2)
{
    return fitClusters(vertex, hit1, hit2);
}
} // namespace LightHelixFit

#endif
