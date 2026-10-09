#include "surface_network.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double smallest_area = 1.0e-10; // m^2, floor that keeps a closed section finite
constexpr double fastest_jet = 60.0;      // m/s

double harmonic_crossing(double first, double second) noexcept {
    double a = std::max(first, smallest_area);
    double b = std::max(second, smallest_area);
    double crossing = 2.0 / (1.0 / a + 1.0 / b);
    return crossing;
}

}

CavityGeometry single_tube(std::span<const double> station_area, double length, std::size_t ring_pieces) {
    CavityGeometry geometry{};
    TubeGeometry airway{};
    airway.stations = station_area.size();
    airway.length = length;
    airway.ring_pieces = ring_pieces;
    airway.end = TubeEnd::radiating;
    geometry.tubes.push_back(airway);
    geometry.station_area.assign(station_area.begin(), station_area.end());
    return geometry;
}

void SurfaceNetwork::build(const CavityGeometry& geometry, double sample_rate, const AirConstants& air) {
    if (geometry.tubes.empty() || geometry.tubes[0].stations < 2 || geometry.tubes[0].parent != no_parent) {
        throw std::invalid_argument("SurfaceNetwork needs an airway tube of at least 2 stations first");
    }
    air_ = air;
    sample_rate_ = sample_rate;

    // Stations of every tube, in tube order.
    stations_.clear();
    tube_first_station_.clear();
    for (std::size_t tube = 0; tube < geometry.tubes.size(); ++tube) {
        const TubeGeometry& shape = geometry.tubes[tube];
        if (shape.stations < 1 || shape.ring_pieces < 1) {
            throw std::invalid_argument("SurfaceNetwork tube needs at least one station and one piece per ring");
        }
        if (tube > 0 && (shape.parent >= tube || shape.parent_station >= geometry.tubes[shape.parent].stations)) {
            throw std::invalid_argument("SurfaceNetwork tube must open from a station of an earlier tube");
        }
        tube_first_station_.push_back(stations_.size());
        for (std::size_t index = 0; index < shape.stations; ++index) {
            Station station{};
            station.tube = tube;
            station.pieces = shape.ring_pieces;
            station.spacing = shape.length / static_cast<double>(shape.stations);
            station.perimeter_factor = shape.perimeter_factor;
            station.yielding_fraction = shape.yielding_fraction;
            station.solid = shape.solid;
            station.cross_paths = shape.solid && shape.cross_paths;
            station.axial_share = station.cross_paths ? std::clamp(shape.axial_share, 0.05, 1.0) : 1.0;
            station.edge_tension = shape.solid ? shape.edge_tension : 0.0;
            station.cross_timing = shape.cross_timing;
            station.fixed_width = shape.width;
            if (shape.solid) {
                double along = (static_cast<double>(index) + 0.5) * station.spacing;
                place_station(shape, along, station);
            }
            station.parallel = shape.parallel;
            if (index > 0) {
                station.upstream = stations_.size() - 1;
            } else if (tube > 0 && !shape.parallel) {
                station.upstream = tube_first_station_[shape.parent] + shape.parent_station;
            }
            stations_.push_back(station);
        }
        for (std::size_t index = 0; index + 1 < shape.stations; ++index) {
            stations_[tube_first_station_[tube] + index].downstream = tube_first_station_[tube] + index + 1;
        }
    }
    beside_.assign(stations_.size(), no_parent);
    for (std::size_t tube = 0; tube < geometry.tubes.size(); ++tube) {
        const TubeGeometry& shape = geometry.tubes[tube];
        if (!shape.parallel) {
            continue;
        }
        if (shape.parent_station + shape.stations > geometry.tubes[shape.parent].stations) {
            throw std::invalid_argument("SurfaceNetwork parallel tube must lie beside stations of its parent");
        }
        for (std::size_t index = 0; index < shape.stations; ++index) {
            beside_[tube_first_station_[tube] + index] = tube_first_station_[shape.parent] + shape.parent_station + index;
        }
    }
    opening_target_.assign(stations_.size(), 0.0);
    opening_now_.assign(stations_.size(), 0.0);
    opening_step_.assign(stations_.size(), 0.0);
    opening_first_.assign(stations_.size(), 0.0);
    opening_area_.assign(stations_.size(), 0.0);
    if (geometry.station_area.size() != stations_.size()) {
        throw std::invalid_argument("SurfaceNetwork station_area must cover every station of every tube");
    }

    // Pieces: the glottis, one ring per station, one aperture per radiating tube.
    pieces_.clear();
    piece_cos_.assign(1, 1.0);
    piece_sin_.assign(1, 0.0);
    SurfacePiece glottis{};
    glottis.normal = Vec3{.x = 1.0, .y = 0.0, .z = 0.0};
    glottis.role = SurfaceRole::glottis;
    glottis.station = 0;
    glottis.yielding = false;
    pieces_.push_back(glottis);
    glottis_ = 0;
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        Station& station = stations_[flat];
        const TubeGeometry& shape = geometry.tubes[station.tube];
        std::size_t index = flat - tube_first_station_[station.tube];
        double x = (static_cast<double>(index) + 0.5) * station.spacing;
        double radius = std::sqrt(std::max(geometry.station_area[flat], smallest_area) / pi);
        double stagger = (index % 2 == 1) ? pi / static_cast<double>(station.pieces) : 0.0;
        if (station.solid) {
            stagger = pi / static_cast<double>(station.pieces);
        }
        station.first_piece = pieces_.size();
        for (std::size_t piece = 0; piece < station.pieces; ++piece) {
            double angle = 2.0 * pi * static_cast<double>(piece) / static_cast<double>(station.pieces) + stagger;
            piece_cos_.push_back(std::cos(angle));
            piece_sin_.push_back(std::sin(angle));
            SurfacePiece wall{};
            wall.position = Vec3{.x = x, .y = radius * std::cos(angle), .z = radius * std::sin(angle)};
            wall.normal = Vec3{.x = 0.0, .y = -std::cos(angle), .z = -std::sin(angle)};
            wall.role = SurfaceRole::wall;
            wall.station = flat;
            wall.tube = station.tube;
            wall.material = shape.material;
            wall.yielding = shape.yielding;
            pieces_.push_back(wall);
        }
    }
    openings_.clear();
    for (std::size_t tube = 0; tube < geometry.tubes.size(); ++tube) {
        if (geometry.tubes[tube].end != TubeEnd::radiating) {
            continue;
        }
        std::size_t last = tube_first_station_[tube] + geometry.tubes[tube].stations - 1;
        SurfacePiece aperture{};
        aperture.position = Vec3{.x = geometry.tubes[tube].length, .y = 0.0, .z = 0.0};
        aperture.normal = Vec3{.x = -1.0, .y = 0.0, .z = 0.0};
        aperture.role = SurfaceRole::aperture;
        aperture.station = last;
        aperture.tube = tube;
        aperture.yielding = false;
        Opening opening{};
        opening.piece = pieces_.size();
        opening.station = last;
        openings_.push_back(opening);
        pieces_.push_back(aperture);
        piece_cos_.push_back(1.0);
        piece_sin_.push_back(0.0);
    }
    // A closed tube ends half a section beyond its last ring: a rigid end
    // piece there reflects every fluctuation with the tube's full length.
    // Without it a side branch was half a section short (MEASURED: the
    // 2 cm piriform branch's antiresonance at 5.8 kHz instead of near
    // 4.3 kHz, a notch fixed in every vowel).
    std::vector<std::size_t> closed_end_piece(geometry.tubes.size(), no_parent);
    for (std::size_t tube = 0; tube < geometry.tubes.size(); ++tube) {
        if (geometry.tubes[tube].end != TubeEnd::closed) {
            continue;
        }
        std::size_t last = tube_first_station_[tube] + geometry.tubes[tube].stations - 1;
        SurfacePiece end{};
        end.role = SurfaceRole::wall;
        end.station = last;
        end.tube = tube;
        end.yielding = false;
        closed_end_piece[tube] = pieces_.size();
        pieces_.push_back(end);
        piece_cos_.push_back(1.0);
        piece_sin_.push_back(0.0);
    }
    const std::size_t count = pieces_.size();
    opening_of_piece_.assign(count, no_parent);
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        opening_of_piece_[openings_[index].piece] = index;
    }

    // Air paths between neighbouring rings of a tube, between a branch's
    // first ring and the ring it opens from, and from the end pieces into
    // their rings. Weight: the cross-section the path passes through
    // (harmonic mean of the two sections), shared over the ring-to-ring
    // paths; an end piece stands for a whole ring. Delay: the axial advance
    // of the wavefront, at least 1.5 samples so that the integer part reads
    // only values already written and the all-pass fraction stays in
    // [0.5, 1.5).
    links_.clear();
    link_length_.clear();
    link_station_.clear();
    double longest = 2.0;
    const double cross_longest = 60.0; // samples: a chord of 22 cm at 192 kHz with the mesh scaling
    axial_scale_.clear();
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        const Station& station = stations_[flat];
        double share_root = std::sqrt(station.axial_share);
        if (station.parallel) {
            continue; // its paths are normalized (below)
        }
        if (station.upstream != no_parent) {
            const Station& upstream = stations_[station.upstream];
            double advance = 0.5 * (station.spacing + upstream.spacing);
            if (station.solid && upstream.solid && station.tube == upstream.tube) {
                add_paired_links(upstream.first_piece, station.first_piece, station.pieces, advance, share_root, flat,
                                 longest);
            } else {
                add_links(upstream.first_piece, upstream.pieces, station.first_piece, station.pieces, advance,
                          share_root, flat, longest);
            }
        } else if (station.tube == 0) {
            add_links(glottis_, 1, station.first_piece, station.pieces, 0.5 * station.spacing, share_root, flat,
                      longest);
        }
    }
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        const Station& station = stations_[openings_[index].station];
        add_links(station.first_piece, station.pieces, openings_[index].piece, 1, 0.5 * station.spacing,
                  std::sqrt(station.axial_share), openings_[index].station, longest);
    }
    for (std::size_t tube = 0; tube < geometry.tubes.size(); ++tube) {
        if (closed_end_piece[tube] == no_parent || geometry.tubes[tube].parallel) {
            continue;
        }
        std::size_t last = tube_first_station_[tube] + geometry.tubes[tube].stations - 1;
        const Station& station = stations_[last];
        add_links(station.first_piece, station.pieces, closed_end_piece[tube], 1, 0.5 * station.spacing,
                  std::sqrt(station.axial_share), last, longest);
    }
    // Cross paths: every pair of pieces of a ring, both directions.
    first_cross_ = links_.size();
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        const Station& station = stations_[flat];
        if (!station.cross_paths) {
            continue;
        }
        for (std::size_t first = 0; first < station.pieces; ++first) {
            for (std::size_t second = first + 1; second < station.pieces; ++second) {
                AirLink forward{.from = station.first_piece + first, .to = station.first_piece + second,
                                .whole_delay = 0, .allpass = 0.0, .admittance = 0.0, .scale = 0.0};
                AirLink backward = forward;
                backward.from = forward.to;
                backward.to = forward.from;
                links_.push_back(forward);
                links_.push_back(backward);
                link_length_.push_back(0.0);
                link_length_.push_back(0.0);
                link_station_.push_back(flat);
                link_station_.push_back(flat);
                axial_scale_.push_back(1.0);
                axial_scale_.push_back(1.0);
            }
        }
    }
    // Parallel tubes stand empty much of the time (the space under the
    // tongue exists only while the blade is lifted): their own paths along
    // the tube and into their closed end are normalized as well, or the
    // pressure waves stored in an empty path grow until it fills.
    first_coupling_ = links_.size();
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        const Station& station = stations_[flat];
        if (!station.parallel) {
            continue;
        }
        std::size_t from = no_parent;
        double advance = 0.0;
        if (station.upstream != no_parent) {
            from = stations_[station.upstream].first_piece;
            advance = 0.5 * (station.spacing + stations_[station.upstream].spacing);
        }
        std::size_t end_piece = closed_end_piece[station.tube];
        for (std::size_t pass = 0; pass < 2; ++pass) {
            std::size_t first = from;
            std::size_t second = station.first_piece;
            double length = advance;
            if (pass == 1) {
                bool last = station.downstream == no_parent;
                if (!last || end_piece == no_parent) {
                    continue;
                }
                first = station.first_piece;
                second = end_piece;
                length = 0.5 * station.spacing;
            }
            if (first == no_parent) {
                continue;
            }
            AirLink forward{.from = first, .to = second, .whole_delay = 0, .allpass = 0.0, .admittance = 0.0,
                            .scale = 0.0};
            AirLink backward = forward;
            backward.from = forward.to;
            backward.to = forward.from;
            links_.push_back(forward);
            links_.push_back(backward);
            link_length_.push_back(length);
            link_length_.push_back(length);
            link_station_.push_back(flat);
            link_station_.push_back(flat);
            axial_scale_.push_back(1.0);
            axial_scale_.push_back(1.0);
        }
    }
    first_opening_link_ = links_.size();
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        if (beside_[flat] == no_parent) {
            continue;
        }
        const Station& station = stations_[flat];
        const Station& parent = stations_[beside_[flat]];
        for (std::size_t own = 0; own < station.pieces; ++own) {
            for (std::size_t other = 0; other < parent.pieces; ++other) {
                AirLink forward{.from = station.first_piece + own, .to = parent.first_piece + other,
                                .whole_delay = 0, .allpass = 0.0, .admittance = 0.0, .scale = 0.0};
                AirLink backward = forward;
                backward.from = forward.to;
                backward.to = forward.from;
                links_.push_back(forward);
                links_.push_back(backward);
                link_length_.push_back(0.0);
                link_length_.push_back(0.0);
                link_station_.push_back(flat);
                link_station_.push_back(flat);
                axial_scale_.push_back(1.0);
                axial_scale_.push_back(1.0);
            }
        }
    }
    const std::size_t crossing_links = links_.size() - first_cross_;
    cross_delay_.assign(crossing_links, 1.0);
    cross_delay_step_.assign(crossing_links, 0.0);
    cross_admittance_.assign(crossing_links, 0.0);
    cross_admittance_step_.assign(crossing_links, 0.0);
    root_admittance_.assign(crossing_links, 0.0);
    if (crossing_links > 0) {
        longest = std::max(longest, cross_longest);
    }
    ring_weight_.assign(stations_.size(), 0.0);
    allpass_step_.assign(links_.size(), 0.0);

    incoming_start_.assign(count + 1, 0);
    outgoing_start_.assign(count + 1, 0);
    for (std::size_t link = 0; link < links_.size(); ++link) {
        incoming_start_[links_[link].to + 1] += 1;
        outgoing_start_[links_[link].from + 1] += 1;
    }
    for (std::size_t piece = 0; piece < count; ++piece) {
        incoming_start_[piece + 1] += incoming_start_[piece];
        outgoing_start_[piece + 1] += outgoing_start_[piece];
    }
    incoming_links_.assign(links_.size(), 0);
    outgoing_links_.assign(links_.size(), 0);
    std::vector<std::size_t> incoming_fill(incoming_start_.begin(), incoming_start_.end() - 1);
    std::vector<std::size_t> outgoing_fill(outgoing_start_.begin(), outgoing_start_.end() - 1);
    for (std::size_t link = 0; link < links_.size(); ++link) {
        incoming_links_[incoming_fill[links_[link].to]] = link;
        incoming_fill[links_[link].to] += 1;
        outgoing_links_[outgoing_fill[links_[link].from]] = link;
        outgoing_fill[links_[link].from] += 1;
    }

    wall_admittance_.assign(count, 0.0);
    wall_filter_.assign(count, std::array<double, 5>{});
    wall_state_.assign(count, std::array<double, 4>{});
    wall_history_.assign(count, 0.0);
    for (std::size_t piece = 0; piece < count; ++piece) {
        if (pieces_[piece].role == SurfaceRole::wall && pieces_[piece].yielding) {
            write_wall_filter(piece);
        }
    }

    admittance_sum_.assign(count, 0.0);
    ring_size_ = std::bit_ceil(static_cast<std::size_t>(std::ceil(longest)) + 4);
    ring_mask_ = ring_size_ - 1;
    for (std::size_t link = 0; link < links_.size(); ++link) {
        links_[link].buffer_offset = link * ring_size_;
    }
    history_.assign(links_.size() * ring_size_, 0.0);
    arriving_.assign(links_.size(), 0.0);
    allpass_input_.assign(links_.size(), 0.0);
    allpass_output_.assign(links_.size(), 0.0);
    loss_filter_.assign(links_.size(), std::array<double, 3>{1.0, 0.0, 0.0});
    loss_state_.assign(links_.size(), std::array<double, 2>{});
    glide_left_ = 0;
    weighted_.assign(count, 0.0);
    pressure_.assign(count, 0.0);
    loss_step_.assign(links_.size(), std::array<double, 3>{});
    const std::size_t boundaries = stations_.size() + openings_.size();
    link_boundary_.assign(links_.size(), no_parent);
    piece_boundary_.assign(count, no_parent);
    // Kinetic losses act on the pressure-wave paths only (not on the
    // normalized ones from first_cross_ on).
    for (std::size_t link = 0; link < first_cross_; ++link) {
        const SurfacePiece& from = pieces_[links_[link].from];
        const SurfacePiece& to = pieces_[links_[link].to];
        if (from.role != SurfaceRole::wall) {
            continue;
        }
        std::size_t boundary = no_parent;
        if (to.role == SurfaceRole::wall && stations_[to.station].upstream == from.station) {
            boundary = to.station;
        } else if (to.role == SurfaceRole::aperture && from.station == to.station) {
            boundary = stations_.size() + opening_of_piece_[links_[link].to];
        }
        link_boundary_[link] = boundary;
        if (boundary != no_parent) {
            piece_boundary_[links_[link].to] = boundary;
        }
    }
    boundary_resistance_.assign(boundaries, 0.0);
    boundary_flow_.assign(boundaries, 0.0);
    boundary_smoothed_.assign(boundaries, 0.0);
    smoothing_ = 1.0 - std::exp(-2.0 * pi * 1000.0 / sample_rate);
    boundary_admittance_.assign(boundaries, 0.0);
    boundary_upstream_.assign(boundaries, no_parent);
    boundary_downstream_.assign(boundaries, no_parent);
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        if (stations_[flat].upstream != no_parent) {
            boundary_upstream_[flat] = stations_[flat].upstream;
            boundary_downstream_[flat] = flat;
        }
    }
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        boundary_upstream_[stations_.size() + index] = openings_[index].station;
    }
    station_area_now_.assign(stations_.size(), 0.0);
    station_area_step_.assign(stations_.size(), 0.0);
    geometry_first_.assign(stations_.size(), 0.0);
    geometry_area_.assign(stations_.size(), 0.0);
    motion_smoothing_ = 1.0 - std::exp(-2.0 * pi * 600.0 / sample_rate);
    contact_smoothing_ = 1.0 - std::exp(-2.0 * pi * 150.0 / sample_rate);
    resisted_sum_.assign(count, 0.0);
    turbulence_boundary_ = no_parent;

    // Edges of solid tubes with tension: ring neighbours and axial neighbours.
    piece_wetted_.assign(count, 0.0);
    displacement_.assign(count, 0.0);
    edge_pressure_.assign(count, 0.0);
    injected_.assign(count, 0.0);
    edge_start_.assign(count + 1, 0);
    edge_piece_.clear();
    for (std::size_t piece = 0; piece < count; ++piece) {
        edge_start_[piece] = edge_piece_.size();
        const SurfacePiece& surface = pieces_[piece];
        if (surface.role != SurfaceRole::wall) {
            continue;
        }
        const Station& station = stations_[surface.station];
        if (!station.solid || station.edge_tension <= 0.0 || !surface.yielding) {
            continue;
        }
        std::size_t index = piece - station.first_piece;
        std::size_t previous = station.first_piece + (index + station.pieces - 1) % station.pieces;
        std::size_t next = station.first_piece + (index + 1) % station.pieces;
        if (station.pieces > 2) {
            edge_piece_.push_back(previous);
            edge_piece_.push_back(next);
        } else if (station.pieces == 2) {
            edge_piece_.push_back(next);
        }
        if (station.upstream != no_parent && stations_[station.upstream].tube == surface.tube) {
            edge_piece_.push_back(stations_[station.upstream].first_piece + index);
        }
        if (station.downstream != no_parent) {
            edge_piece_.push_back(stations_[station.downstream].first_piece + index);
        }
    }
    edge_start_[count] = edge_piece_.size();
    edge_gain_.assign(edge_piece_.size(), 0.0);

    reshape(std::span<const double>(geometry.station_area), 1);
    reset();
}

void SurfaceNetwork::add_links(std::size_t first_piece, std::size_t first_count, std::size_t second_piece,
                               std::size_t second_count, double advance, double share_root,
                               std::size_t lining_station, double& longest) {
    const double characteristic = air_.density * air_.sound_speed;
    double delay = std::max(advance, 1.0e-5) / air_.sound_speed * sample_rate_ * share_root;
    delay = std::max(delay, 1.5);
    std::size_t whole = static_cast<std::size_t>(std::floor(delay - 1.5));
    double fraction = delay - 1.0 - static_cast<double>(whole);
    double allpass = (1.0 - fraction) / (1.0 + fraction);
    double scale = share_root / (static_cast<double>(first_count * second_count) * characteristic);
    for (std::size_t first = 0; first < first_count; ++first) {
        for (std::size_t second = 0; second < second_count; ++second) {
            AirLink forward{.from = first_piece + first, .to = second_piece + second, .whole_delay = whole,
                            .allpass = allpass, .admittance = 0.0, .scale = scale};
            AirLink backward = forward;
            backward.from = forward.to;
            backward.to = forward.from;
            links_.push_back(forward);
            links_.push_back(backward);
            link_length_.push_back(advance);
            link_length_.push_back(advance);
            link_station_.push_back(lining_station);
            link_station_.push_back(lining_station);
            axial_scale_.push_back(share_root);
            axial_scale_.push_back(share_root);
        }
    }
    longest = std::max(longest, delay);
}

// Axial paths of a solid tube: piece k of one ring to piece k of the next,
// each carrying its share of the section. The integer part of the delay is
// fixed here from the advance along the outer wall's midline; the all-pass
// then follows each pair's depth (solid_targets), within
// [whole + 1.05, whole + 3].
void SurfaceNetwork::add_paired_links(std::size_t first_piece, std::size_t second_piece, std::size_t count,
                                      double advance, double share_root, std::size_t lining_station,
                                      double& longest) {
    const double characteristic = air_.density * air_.sound_speed;
    double delay = std::max(advance, 1.0e-5) / air_.sound_speed * sample_rate_ * share_root;
    std::size_t whole = static_cast<std::size_t>(std::max(std::floor(delay - 1.6), 0.0));
    double fraction = std::clamp(delay - 1.0 - static_cast<double>(whole), 0.05, 2.0);
    double allpass = (1.0 - fraction) / (1.0 + fraction);
    double scale = share_root / (static_cast<double>(count) * characteristic);
    for (std::size_t piece = 0; piece < count; ++piece) {
        AirLink forward{.from = first_piece + piece, .to = second_piece + piece, .whole_delay = whole,
                        .allpass = allpass, .admittance = 0.0, .scale = scale};
        AirLink backward = forward;
        backward.from = forward.to;
        backward.to = forward.from;
        links_.push_back(forward);
        links_.push_back(backward);
        link_length_.push_back(advance);
        link_length_.push_back(advance);
        link_station_.push_back(lining_station);
        link_station_.push_back(lining_station);
        axial_scale_.push_back(share_root);
        axial_scale_.push_back(share_root);
    }
    longest = std::max(longest, static_cast<double>(whole) + 3.0);
}

// The outer wall at a distance along a solid tube. Straight: the wall runs
// along x and the section reaches up (+z). Bent (the airway): the posterior
// pharyngeal wall rises from the glottis for 0.37 of the length, turns
// through a quarter circle of outer radius 4 cm (the length is measured
// along a 3 cm midline) and runs forward under the palate to the lips;
// widths follow the anatomy (narrow epilarynx, round pharynx, wide flat
// mouth, wider lips), scaled with the tract length.
void SurfaceNetwork::place_station(const TubeGeometry& shape, double along, Station& station) noexcept {
    station.bend_radius = 0.0;
    station.half_width = shape.width > 0.0 ? 0.5 * shape.width : 0.0;
    if (!shape.bent) {
        station.wall = Vec3{.x = along, .y = 0.0, .z = 0.0};
        station.inward = Vec3{.x = 0.0, .y = 0.0, .z = 1.0};
        return;
    }
    const double scale = shape.length / 0.175;
    const double outer = 0.040 * scale;
    const double midline = 0.030 * scale;
    const double pharynx = 0.37 * shape.length;
    const double bend_end = pharynx + 0.5 * pi * midline;
    if (along <= pharynx) {
        station.wall = Vec3{.x = 0.0, .y = 0.0, .z = along};
        station.inward = Vec3{.x = 1.0, .y = 0.0, .z = 0.0};
    } else if (along < bend_end) {
        double angle = (along - pharynx) / midline;
        station.wall = Vec3{.x = outer - outer * std::cos(angle), .y = 0.0, .z = pharynx + outer * std::sin(angle)};
        station.inward = Vec3{.x = std::cos(angle), .y = 0.0, .z = -std::sin(angle)};
        station.bend_radius = outer;
    } else {
        station.wall = Vec3{.x = outer + (along - bend_end), .y = 0.0, .z = pharynx + outer};
        station.inward = Vec3{.x = 0.0, .y = 0.0, .z = -1.0};
    }
    struct WidthKnot final {
        double along{};
        double width{};
    };
    constexpr WidthKnot knots[] = {
        {0.00, 0.012}, {0.10, 0.014}, {0.16, 0.024}, {0.45, 0.030}, {0.55, 0.034}, {0.88, 0.034}, {1.00, 0.040},
    };
    const std::size_t count = sizeof(knots) / sizeof(knots[0]);
    double fraction = along / shape.length;
    double width = knots[count - 1].width;
    if (fraction <= knots[0].along) {
        width = knots[0].width;
    } else {
        for (std::size_t knot = 1; knot < count; ++knot) {
            if (fraction <= knots[knot].along) {
                double t = (fraction - knots[knot - 1].along) / (knots[knot].along - knots[knot - 1].along);
                double eased = t * t * (3.0 - 2.0 * t);
                width = knots[knot - 1].width + eased * (knots[knot].width - knots[knot - 1].width);
                break;
            }
        }
    }
    station.half_width = 0.5 * width * scale;
}

double SurfaceNetwork::path_factor(const Station& station, double depth, double area, std::size_t flat) const noexcept {
    if (station.bend_radius <= 0.0) {
        return 1.0;
    }
    double a = half_width_of(flat, area);
    double c = std::max(area, smallest_area) / (pi * a);
    double factor = (station.bend_radius - depth) / std::max(station.bend_radius - c, 0.2 * station.bend_radius);
    return factor;
}

double SurfaceNetwork::half_width_of(std::size_t station, double area) const noexcept {
    double half = stations_[station].half_width;
    if (half <= 0.0) {
        half = std::sqrt(std::max(area, smallest_area) / pi);
    }
    return half;
}

// Position of a wall piece of a solid ring for the given section area: on
// the ellipse of half width a and half depth c = A / (pi a) that touches the
// outer wall, at angle phi (0 at the wall). Returns the depth of the piece
// from the outer wall; position and inward normal go to the outputs.
double SurfaceNetwork::piece_place(std::size_t piece, double area, Vec3& position, Vec3& inward) const noexcept {
    const SurfacePiece& surface = pieces_[piece];
    const Station& station = stations_[surface.station];
    double a = half_width_of(surface.station, area);
    double c = std::max(area, smallest_area) / (pi * a);
    double cosine = piece_cos_[piece];
    double sine = piece_sin_[piece];
    double depth = c * (1.0 - cosine);
    double lateral = a * sine;
    position = Vec3{.x = station.wall.x + depth * station.inward.x, .y = lateral,
                    .z = station.wall.z + depth * station.inward.z};
    // Outward normal of the ellipse in (lateral, depth-from-centre) is
    // proportional to (sin/a, -cos/c); the cavity side is its opposite.
    double normal_lateral = -sine / a;
    double normal_depth = cosine / c;
    double norm = std::sqrt(normal_lateral * normal_lateral + normal_depth * normal_depth);
    normal_lateral /= norm;
    normal_depth /= norm;
    inward = Vec3{.x = normal_depth * station.inward.x, .y = normal_lateral, .z = normal_depth * station.inward.z};
    return depth;
}

double SurfaceNetwork::piece_arc(std::size_t piece, double area) const noexcept {
    const SurfacePiece& surface = pieces_[piece];
    const Station& station = stations_[surface.station];
    double a = half_width_of(surface.station, area);
    double c = std::max(area, smallest_area) / (pi * a);
    double cosine = piece_cos_[piece];
    double sine = piece_sin_[piece];
    double speed = std::sqrt(a * a * cosine * cosine + c * c * sine * sine);
    double arc = speed * 2.0 * pi / static_cast<double>(station.pieces);
    return arc;
}

void SurfaceNetwork::write_wall_filter(std::size_t piece) noexcept {
    // z(s) = (m s^2 + r s + k) / (rho c s); R = (z - 1)/(z + 1)
    //   = (m s^2 + (r - rho c) s + k) / (m s^2 + (r + rho c) s + k).
    const WallMaterial& material = pieces_[piece].material;
    double characteristic = air_.density * air_.sound_speed;
    double c = 2.0 * sample_rate_;
    double m = material.mass;
    double k = material.stiffness;
    double low = material.resistance - characteristic;
    double high = material.resistance + characteristic;
    // Bilinear: s = c (1 - z^-1)/(1 + z^-1); multiply through by (1 + z^-1)^2.
    double n0 = m * c * c + low * c + k;
    double n1 = -2.0 * m * c * c + 2.0 * k;
    double n2 = m * c * c - low * c + k;
    double d0 = m * c * c + high * c + k;
    double d1 = -2.0 * m * c * c + 2.0 * k;
    double d2 = m * c * c - high * c + k;
    wall_filter_[piece] = std::array<double, 5>{n0 / d0, n1 / d0, n2 / d0, d1 / d0, d2 / d0};
}

// Radiation port: R(s) = (s tau (r-1) - r) / (s tau (r+1) + r) with
// tau = 8 a / (3 pi c), r = 128 / (9 pi^2), bilinear transform. It follows
// the opening's area.
void SurfaceNetwork::write_radiation(Opening& opening, double area) noexcept {
    double radius = std::sqrt(std::max(area, smallest_area) / pi);
    double tau = 8.0 * radius / (3.0 * pi * air_.sound_speed);
    double resistance = 128.0 / (9.0 * pi * pi);
    // A vanishing opening tends to a pure short (b0 -> -1), where the
    // instantaneous solve divides by 1 + b0. The paths into a nearly sealed
    // opening carry almost no air, so its filter is held at k >= 1 (radius
    // >= 0.34 mm at 192 kHz) without audible effect.
    double k = std::max(2.0 * sample_rate_ * tau, 1.0);
    double denominator = k * (resistance + 1.0) + resistance;
    opening.b0 = (k * (resistance - 1.0) - resistance) / denominator;
    opening.b1 = (-k * (resistance - 1.0) - resistance) / denominator;
    opening.a1 = (-k * (resistance + 1.0) + resistance) / denominator;
}

// Loss over a path of length dx through a section of area A: the
// boundary-layer attenuation alpha(f) dx with
//   alpha(f) = (P / 2A) sqrt(pi f mu / rho) / c (1 + heat factor),
// P / 2A = perimeter_factor / radius. Modelled as the shelf
// H(s) = (1 + s/wz) / (1 + s/wp), wp = 2 pi 3 kHz, whose loss
// eps x^2/(1 + x^2) (x = f/3 kHz) follows sqrt(f) at and above 3 kHz.
std::array<double, 3> SurfaceNetwork::boundary_loss(std::size_t link, double crossing) const noexcept {
    constexpr double knee = 3000.0;
    const Station& lining = stations_[link_station_[link]];
    double radius = std::sqrt(std::max(crossing, smallest_area) / pi);
    double perimeter_over_area = lining.perimeter_factor / radius; // P / 2A
    if (lining.solid) {
        // The ellipse's own perimeter (Ramanujan's approximation).
        double a = half_width_of(link_station_[link], crossing);
        double c = std::max(crossing, smallest_area) / (pi * a);
        double perimeter = pi * (3.0 * (a + c) - std::sqrt((3.0 * a + c) * (a + 3.0 * c)));
        perimeter_over_area = perimeter / (2.0 * std::max(crossing, smallest_area));
    }
    double alpha = perimeter_over_area * std::sqrt(pi * knee * air_.viscosity / air_.density)
        / air_.sound_speed * (1.0 + air_.heat_loss_factor);
    double epsilon = std::min(2.0 * alpha * link_length_[link], 0.5);
    double pole = 2.0 * pi * knee;
    double zero = pole * std::exp(epsilon); // high-frequency gain wp / wz = exp(-eps)
    // Bilinear transform with the knee prewarped.
    double c = pole / std::tan(pole / (2.0 * sample_rate_));
    double n0 = 1.0 + c / zero;
    double n1 = 1.0 - c / zero;
    double d0 = 1.0 + c / pole;
    double d1 = 1.0 - c / pole;
    std::array<double, 3> coefficients{n0 / d0, n1 / d0, d1 / d0};
    return coefficients;
}

double SurfaceNetwork::crossing_area(std::size_t link, std::span<const double> station_area) const noexcept {
    double crossing = harmonic_crossing(station_area[pieces_[links_[link].from].station],
                                        station_area[pieces_[links_[link].to].station]);
    return crossing;
}

void SurfaceNetwork::reshape(std::span<const double> station_area, std::size_t glide_samples) noexcept {
    const std::size_t glide = std::max<std::size_t>(glide_samples, 1);
    // Finish any glide in progress at its target before starting the next.
    while (glide_left_ > 0) {
        advance_glide();
    }
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        double area = std::max(station_area[flat], smallest_area);
        if (glide == 1) {
            station_area_now_[flat] = area;
            station_area_step_[flat] = 0.0;
            geometry_first_[flat] = area;
            geometry_area_[flat] = area;
        } else {
            station_area_step_[flat] = (area - station_area_now_[flat]) / static_cast<double>(glide);
        }
        if (glide == 1) {
            opening_now_[flat] = opening_target_[flat];
            opening_first_[flat] = opening_target_[flat];
            opening_area_[flat] = opening_target_[flat];
            opening_step_[flat] = 0.0;
        } else {
            opening_step_[flat] = (opening_target_[flat] - opening_now_[flat]) / static_cast<double>(glide);
        }
    }
    // Boundary-layer loss filters of the pressure-wave paths (the
    // normalized ones carry none).
    for (std::size_t link = 0; link < first_cross_; ++link) {
        double crossing = crossing_area(link, station_area);
        std::array<double, 3> loss = boundary_loss(link, crossing);
        if (glide == 1) {
            loss_filter_[link] = loss;
            loss_step_[link] = std::array<double, 3>{};
        } else {
            for (std::size_t index = 0; index < 3; ++index) {
                loss_step_[link][index] = (loss[index] - loss_filter_[link][index]) / static_cast<double>(glide);
            }
        }
    }
    solid_targets(station_area, glide);
    glide_left_ = glide == 1 ? 0 : glide;
    if (glide == 1) {
        refresh_admittances();
    }
}

// Targets of the solid tubes for new section areas: the all-pass of every
// axial path from its pieces' depth (in a bend the inner side is shorter),
// the travel time and air of every cross path from its chord, and the pull
// of every edge.
//
// Cross paths share a ring's air by the measure of the rays between their
// pieces (two-dimensional form factor times chord: s_i s_j cos_i cos_j),
// which counts every direction across the section equally; each path holds
// that share of (1 - axial share) of the ring's volume. Its travel time is
// the chord over c, shortened by sqrt((1 - axial share) / 2) as in a
// two-dimensional mesh; its admittance is its volume over rho c^2 per unit
// of travel time.
void SurfaceNetwork::solid_targets(std::span<const double> station_area, std::size_t glide) noexcept {
    const double c = air_.sound_speed;
    const double stiffness = air_.density * c * c;
    const double inverse_glide = 1.0 / static_cast<double>(glide);
    for (std::size_t link = 0; link < first_cross_; ++link) {
        const SurfacePiece& from = pieces_[links_[link].from];
        const SurfacePiece& to = pieces_[links_[link].to];
        if (from.role != SurfaceRole::wall || to.role != SurfaceRole::wall) {
            continue;
        }
        const Station& first = stations_[from.station];
        const Station& second = stations_[to.station];
        if (!first.solid || !second.solid || from.tube != to.tube || first.bend_radius + second.bend_radius <= 0.0) {
            continue;
        }
        Vec3 position{};
        Vec3 inward{};
        double depth_from = piece_place(links_[link].from, station_area[from.station], position, inward);
        double depth_to = piece_place(links_[link].to, station_area[to.station], position, inward);
        double factor = 0.5 * (path_factor(first, depth_from, station_area[from.station], from.station)
                               + path_factor(second, depth_to, station_area[to.station], to.station));
        double advance = 0.5 * (first.spacing + second.spacing) * std::max(factor, 0.2);
        double delay = advance / c * sample_rate_ * axial_scale_[link];
        double fraction = std::clamp(delay - 1.0 - static_cast<double>(links_[link].whole_delay), 0.05, 2.0);
        double target = (1.0 - fraction) / (1.0 + fraction);
        if (glide == 1) {
            links_[link].allpass = target;
            allpass_step_[link] = 0.0;
        } else {
            allpass_step_[link] = (target - links_[link].allpass) * inverse_glide;
        }
    }
    const std::size_t crossing_links = first_coupling_ - first_cross_;
    if (crossing_links > 0) {
        for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
            ring_weight_[flat] = 0.0;
        }
        // First pass: ray measure of every pair (stored as the admittance
        // target for now), summed per ring.
        for (std::size_t index = 0; index < crossing_links; index += 2) {
            std::size_t link = first_cross_ + index;
            std::size_t from = links_[link].from;
            std::size_t to = links_[link].to;
            std::size_t flat = pieces_[from].station;
            double area = std::max(station_area[flat], smallest_area);
            Vec3 first_position{};
            Vec3 first_inward{};
            Vec3 second_position{};
            Vec3 second_inward{};
            piece_place(from, area, first_position, first_inward);
            piece_place(to, area, second_position, second_inward);
            double dx = second_position.x - first_position.x;
            double dy = second_position.y - first_position.y;
            double dz = second_position.z - first_position.z;
            double chord = std::sqrt(dx * dx + dy * dy + dz * dz);
            double weight = 0.0;
            if (chord > 1.0e-9) {
                double first_facing = (dx * first_inward.x + dy * first_inward.y + dz * first_inward.z) / chord;
                double second_facing = -(dx * second_inward.x + dy * second_inward.y + dz * second_inward.z) / chord;
                weight = piece_arc(from, area) * piece_arc(to, area) * std::max(first_facing, 0.0)
                    * std::max(second_facing, 0.0);
            }
            cross_admittance_step_[index] = weight;
            cross_delay_step_[index] = chord;
            ring_weight_[flat] += weight;
        }
        // Second pass: shares, travel times and admittances.
        for (std::size_t index = 0; index < crossing_links; index += 2) {
            std::size_t link = first_cross_ + index;
            std::size_t flat = pieces_[links_[link].from].station;
            const Station& station = stations_[flat];
            double area = std::max(station_area[flat], smallest_area);
            double mesh = std::sqrt(0.5 * (1.0 - station.axial_share)) * station.cross_timing;
            double weight = cross_admittance_step_[index];
            double chord = cross_delay_step_[index];
            double share = ring_weight_[flat] > 0.0 ? weight / ring_weight_[flat] : 0.0;
            double volume = share * (1.0 - station.axial_share) * area * station.spacing;
            double delay = std::clamp(mesh * chord / c * sample_rate_, 1.0, 60.0);
            // Per unit of the ring's area: refresh_admittances multiplies by
            // the smoothed area each sample.
            double admittance = volume / (stiffness * delay / sample_rate_) / area;
            for (std::size_t direction = 0; direction < 2; ++direction) {
                std::size_t slot = index + direction;
                if (glide == 1) {
                    cross_delay_[slot] = delay;
                    cross_admittance_[slot] = admittance;
                    cross_delay_step_[slot] = 0.0;
                    cross_admittance_step_[slot] = 0.0;
                } else {
                    cross_delay_step_[slot] = (delay - cross_delay_[slot]) * inverse_glide;
                    cross_admittance_step_[slot] = (admittance - cross_admittance_[slot]) * inverse_glide;
                }
            }
        }
    }
    // Couplings: the face through which a parallel station opens into the
    // parent station beside it.
    const double characteristic = air_.density * c;
    for (std::size_t link = first_coupling_; link < first_opening_link_; ++link) {
        // Along a parallel tube: the harmonic crossing of its two sections.
        std::size_t slot = link - first_cross_;
        double first_area = station_area[pieces_[links_[link].from].station];
        double second_area = station_area[pieces_[links_[link].to].station];
        double admittance = harmonic_crossing(first_area, second_area) / characteristic;
        double delay = std::clamp(link_length_[link] / c * sample_rate_, 1.0, 60.0);
        if (glide == 1) {
            cross_delay_[slot] = delay;
            cross_admittance_[slot] = admittance;
            cross_delay_step_[slot] = 0.0;
            cross_admittance_step_[slot] = 0.0;
        } else {
            cross_delay_step_[slot] = (delay - cross_delay_[slot]) * inverse_glide;
            cross_admittance_step_[slot] = (admittance - cross_admittance_[slot]) * inverse_glide;
        }
    }
    for (std::size_t link = first_opening_link_; link < links_.size(); ++link) {
        std::size_t slot = link - first_cross_;
        std::size_t own = pieces_[links_[link].from].station;
        std::size_t other = pieces_[links_[link].to].station;
        std::size_t flat = beside_[own] != no_parent ? own : other;
        std::size_t parent = beside_[flat];
        double share = 1.0 / static_cast<double>(stations_[parent].pieces * stations_[flat].pieces);
        double admittance = opening_target_[flat] / characteristic * share;
        // The face opens straight into the section beside it: half a
        // station of travel.
        double delay = std::clamp(0.5 * stations_[flat].spacing / c * sample_rate_, 1.0, 60.0);
        if (glide == 1) {
            cross_delay_[slot] = delay;
            cross_admittance_[slot] = admittance;
            cross_delay_step_[slot] = 0.0;
            cross_admittance_step_[slot] = 0.0;
        } else {
            cross_delay_step_[slot] = (delay - cross_delay_[slot]) * inverse_glide;
            cross_admittance_step_[slot] = (admittance - cross_admittance_[slot]) * inverse_glide;
        }
    }
    // Edges: e_i = sum T g (x_j - x_i) / S_i, g = shared edge length over
    // the distance between the pieces' centres.
    for (std::size_t piece = 0; piece < pieces_.size(); ++piece) {
        if (edge_start_[piece] == edge_start_[piece + 1]) {
            continue;
        }
        const Station& station = stations_[pieces_[piece].station];
        double own_arc = piece_arc(piece, station_area[pieces_[piece].station]);
        double surface = std::max(own_arc * station.spacing, 1.0e-6);
        for (std::size_t position = edge_start_[piece]; position < edge_start_[piece + 1]; ++position) {
            std::size_t other = edge_piece_[position];
            const Station& other_station = stations_[pieces_[other].station];
            double other_arc = piece_arc(other, station_area[pieces_[other].station]);
            double geometry = 0.0;
            if (pieces_[other].station == pieces_[piece].station) {
                geometry = station.spacing / std::max(0.5 * (own_arc + other_arc), 1.0e-4);
            } else {
                geometry = 0.5 * (own_arc + other_arc) / (0.5 * (station.spacing + other_station.spacing));
            }
            edge_gain_[position] = station.edge_tension * geometry / surface;
        }
    }
}

// Every admittance from the current geometry areas: walls from the wetted
// perimeter, paths from the harmonic crossing, each radiation port matched to
// the paths arriving at its aperture.
void SurfaceNetwork::refresh_admittances() noexcept {
    const std::size_t count = pieces_.size();
    const double characteristic = air_.density * air_.sound_speed;
    for (std::size_t piece = 0; piece < count; ++piece) {
        double wall = 0.0;
        const SurfacePiece& surface = pieces_[piece];
        if (surface.role == SurfaceRole::wall && surface.yielding) {
            const Station& station = stations_[surface.station];
            double wetted = 0.0;
            if (station.solid) {
                wetted = station.yielding_fraction * piece_arc(piece, geometry_area_[surface.station])
                    * station.spacing;
            } else {
                double radius = std::sqrt(geometry_area_[surface.station] / pi);
                wetted = station.perimeter_factor * station.yielding_fraction * 2.0 * pi * radius
                    * station.spacing / static_cast<double>(station.pieces);
            }
            piece_wetted_[piece] = wetted;
            wall = wetted / characteristic;
        }
        wall_admittance_[piece] = wall;
        admittance_sum_[piece] = wall;
    }
    for (std::size_t link = 0; link < first_cross_; ++link) {
        double crossing = harmonic_crossing(geometry_area_[pieces_[links_[link].from].station],
                                            geometry_area_[pieces_[links_[link].to].station]);
        double admittance = crossing * links_[link].scale;
        links_[link].admittance = admittance;
        admittance_sum_[links_[link].to] += admittance;
    }
    // Normalized paths follow the smoothed geometry sample by sample, as the
    // pressure-wave paths do: gliding linearly per block, their admittances
    // changed slope at every block edge, and a fast-moving tongue tip
    // clicked at each (MEASURED: tip onsets 7-13x louder with the space
    // under the tongue than without it).
    // Cross paths: admittance per unit of their ring's area (glided per
    // block) times the ring's smoothed area.
    for (std::size_t link = first_cross_; link < first_coupling_; ++link) {
        std::size_t slot = link - first_cross_;
        double admittance = std::max(cross_admittance_[slot], 0.0) * geometry_area_[pieces_[links_[link].from].station];
        links_[link].admittance = admittance;
        root_admittance_[slot] = std::sqrt(admittance);
        admittance_sum_[links_[link].to] += admittance;
    }
    // Along a parallel tube: the harmonic crossing of its two sections.
    for (std::size_t link = first_coupling_; link < first_opening_link_; ++link) {
        std::size_t slot = link - first_cross_;
        double admittance = harmonic_crossing(geometry_area_[pieces_[links_[link].from].station],
                                              geometry_area_[pieces_[links_[link].to].station])
            / characteristic;
        links_[link].admittance = admittance;
        root_admittance_[slot] = std::sqrt(admittance);
        admittance_sum_[links_[link].to] += admittance;
    }
    // Openings of parallel stations into the stations beside them, shared
    // over the pieces of both rings.
    for (std::size_t link = first_opening_link_; link < links_.size(); ++link) {
        std::size_t slot = link - first_cross_;
        std::size_t own = pieces_[links_[link].from].station;
        std::size_t other = pieces_[links_[link].to].station;
        std::size_t flat = beside_[own] != no_parent ? own : other;
        double share = 1.0 / static_cast<double>(stations_[beside_[flat]].pieces * stations_[flat].pieces);
        double admittance = opening_area_[flat] / characteristic * share;
        links_[link].admittance = admittance;
        root_admittance_[slot] = std::sqrt(admittance);
        admittance_sum_[links_[link].to] += admittance;
    }
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        Opening& opening = openings_[index];
        opening.admittance = admittance_sum_[opening.piece];
        admittance_sum_[opening.piece] += opening.admittance;
        opening.area = geometry_area_[opening.station];
        write_radiation(opening, opening.area);
    }
}

// One sample of motion. The station areas glide linearly to the block's
// targets and then pass a two-pole smoothing (about 0.5 ms), so the cavity
// moves with continuous velocity: with linear glides alone the rate of
// closure jumped at every block edge, and the microphone, which hears the
// rate of change of flow, heard each jump as a click (MEASURED in whisper,
// 30x the transients of a human whisper).
void SurfaceNetwork::advance_glide() noexcept {
    if (glide_left_ > 0) {
        for (std::size_t link = 0; link < first_cross_; ++link) {
            links_[link].allpass += allpass_step_[link];
        }
        const std::size_t crossing_links = links_.size() - first_cross_;
        for (std::size_t index = 0; index < crossing_links; ++index) {
            cross_delay_[index] += cross_delay_step_[index];
            cross_admittance_[index] += cross_admittance_step_[index];
        }
        for (std::size_t link = 0; link < first_cross_; ++link) {
            std::array<double, 3>& loss = loss_filter_[link];
            const std::array<double, 3>& step = loss_step_[link];
            loss[0] += step[0];
            loss[1] += step[1];
            loss[2] += step[2];
        }
        for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
            station_area_now_[flat] += station_area_step_[flat];
            opening_now_[flat] += opening_step_[flat];
        }
        glide_left_ -= 1;
    }
    // Soft tissue meets gradually: the last fraction of a millimetre of a
    // closing (or opening) gap moves slowly. Below 0.2 cm^2 the smoothing
    // slows from 600 Hz toward 150 Hz at contact; sealed at full speed, a
    // closure cut a whisper's airflow off at once (a thump 3x the vowel).
    const double contact_area = 0.2e-4;
    for (std::size_t flat = 0; flat < stations_.size(); ++flat) {
        double nearness = std::clamp(1.0 - geometry_area_[flat] / contact_area, 0.0, 1.0);
        double smoothing = motion_smoothing_ + nearness * (contact_smoothing_ - motion_smoothing_);
        geometry_first_[flat] += smoothing * (station_area_now_[flat] - geometry_first_[flat]);
        geometry_area_[flat] += smoothing * (geometry_first_[flat] - geometry_area_[flat]);
        geometry_area_[flat] = std::max(geometry_area_[flat], smallest_area);
        opening_first_[flat] += motion_smoothing_ * (opening_now_[flat] - opening_first_[flat]);
        opening_area_[flat] += motion_smoothing_ * (opening_first_[flat] - opening_area_[flat]);
    }
    refresh_admittances();
}

std::size_t SurfaceNetwork::exit_boundary(std::size_t station) const noexcept {
    if (station >= stations_.size()) {
        return no_parent;
    }
    if (stations_[station].downstream != no_parent) {
        return stations_[station].downstream;
    }
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        if (openings_[index].station == station) {
            return stations_.size() + index;
        }
    }
    return no_parent;
}

void SurfaceNetwork::set_turbulence_exit(std::size_t station, std::size_t fade_samples) noexcept {
    std::size_t boundary = exit_boundary(station);
    if (boundary == turbulence_boundary_) {
        return;
    }
    fading_boundary_ = turbulence_weight_ >= 0.5 ? turbulence_boundary_ : fading_boundary_;
    turbulence_boundary_ = boundary;
    turbulence_weight_ = 0.0;
    turbulence_fade_ = 1.0 / static_cast<double>(std::max<std::size_t>(fade_samples, 1));
    if (fading_boundary_ == boundary) {
        fading_boundary_ = no_parent;
        turbulence_weight_ = 1.0;
    }
}

double SurfaceNetwork::exit_flow(std::size_t station) const noexcept {
    std::size_t boundary = exit_boundary(station);
    if (boundary == no_parent) {
        return 0.0;
    }
    return boundary_flow_[boundary];
}

void SurfaceNetwork::reset() noexcept {
    for (std::size_t index = 0; index < history_.size(); ++index) {
        history_[index] = 0.0;
    }
    for (std::size_t link = 0; link < links_.size(); ++link) {
        arriving_[link] = 0.0;
        allpass_input_[link] = 0.0;
        allpass_output_[link] = 0.0;
        loss_state_[link] = std::array<double, 2>{};
    }
    for (std::size_t piece = 0; piece < pieces_.size(); ++piece) {
        pressure_[piece] = 0.0;
        weighted_[piece] = 0.0;
        displacement_[piece] = 0.0;
        edge_pressure_[piece] = 0.0;
        wall_state_[piece] = std::array<double, 4>{};
        wall_history_[piece] = 0.0;
    }
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        Opening& opening = openings_[index];
        opening.previous_into = 0.0;
        opening.previous_out = 0.0;
        opening.history = 0.0;
        opening.flow = 0.0;
        opening.previous_flow = 0.0;
    }
    for (std::size_t boundary = 0; boundary < boundary_flow_.size(); ++boundary) {
        boundary_flow_[boundary] = 0.0;
        boundary_smoothed_[boundary] = 0.0;
        boundary_resistance_[boundary] = 0.0;
    }
    write_index_ = 0;
    wall_flow_ = 0.0;
}

void SurfaceNetwork::gather() noexcept {
    // The geometry moves first, so this sample's arrivals, junctions and
    // emissions all see the same admittances.
    advance_glide();
    // Every arriving fluctuation: the value emitted whole_delay + 1 samples
    // ago, through the path's fractional all-pass and its boundary-layer loss.
    for (std::size_t link = 0; link < first_cross_; ++link) {
        const AirLink& path = links_[link];
        double delayed = history_[path.buffer_offset + ((write_index_ - path.whole_delay - 1) & ring_mask_)];
        double value = path.allpass * delayed + allpass_input_[link] - path.allpass * allpass_output_[link];
        allpass_input_[link] = delayed;
        allpass_output_[link] = value;
        const std::array<double, 3>& loss = loss_filter_[link];
        std::array<double, 2>& state = loss_state_[link];
        double lossy = loss[0] * value + loss[1] * state[0] - loss[2] * state[1];
        state[0] = value;
        state[1] = lossy;
        arriving_[link] = lossy;
    }
    // Cross paths and openings carry power-normalized waves (sqrt(Y) times
    // the pressure wave): their admittance swings to zero and back as
    // sections close and openings move, and a pressure wave stored in a
    // path that no longer acts on its ends accumulates without bound (each
    // round trip adds the pressure difference of its ends; MEASURED: the
    // first opening under the tongue blew the voice up). A normalized wave
    // in a path without air only circulates, and a changing sqrt(Y) neither
    // creates nor destroys its energy.
    // The value emitted delay samples ago, linearly interpolated.
    for (std::size_t link = first_cross_; link < links_.size(); ++link) {
        double delay = cross_delay_[link - first_cross_];
        double whole = std::floor(delay);
        double fraction = delay - whole;
        std::size_t back = static_cast<std::size_t>(whole);
        std::size_t offset = links_[link].buffer_offset;
        double nearer = history_[offset + ((write_index_ - back) & ring_mask_)];
        double farther = history_[offset + ((write_index_ - back - 1) & ring_mask_)];
        arriving_[link] = nearer + fraction * (farther - nearer);
    }
    // The pull of the edges on each wall piece, from the last displacements.
    const std::size_t count = pieces_.size();
    for (std::size_t piece = 0; piece < count; ++piece) {
        double pull = 0.0;
        double own = displacement_[piece];
        for (std::size_t position = edge_start_[piece]; position < edge_start_[piece + 1]; ++position) {
            pull += edge_gain_[position] * (displacement_[edge_piece_[position]] - own);
        }
        edge_pressure_[piece] = pull;
    }
    for (std::size_t piece = 0; piece < count; ++piece) {
        double weighted = 0.0;
        double resisted = 0.0;
        std::size_t boundary = piece_boundary_[piece];
        if (boundary != no_parent) {
            // A resisted path acts as admittance Y / (1 + R Y).
            double through = 1.0 / (1.0 + boundary_resistance_[boundary]);
            for (std::size_t position = incoming_start_[piece]; position < incoming_start_[piece + 1]; ++position) {
                std::size_t link = incoming_links_[position];
                double contribution = link < first_cross_ ? links_[link].admittance * arriving_[link]
                                                          : root_admittance_[link - first_cross_] * arriving_[link];
                if (link_boundary_[link] != no_parent) {
                    contribution *= through;
                    resisted += links_[link].admittance;
                }
                weighted += contribution;
            }
        } else {
            for (std::size_t position = incoming_start_[piece]; position < incoming_start_[piece + 1]; ++position) {
                std::size_t link = incoming_links_[position];
                if (link < first_cross_) {
                    weighted += links_[link].admittance * arriving_[link];
                } else {
                    weighted += root_admittance_[link - first_cross_] * arriving_[link];
                }
            }
        }
        resisted_sum_[piece] = resisted;
        weighted_[piece] = weighted;
        if (wall_admittance_[piece] > 0.0) {
            const std::array<double, 5>& f = wall_filter_[piece];
            const std::array<double, 4>& state = wall_state_[piece];
            wall_history_[piece] = f[1] * state[0] + f[2] * state[1] - f[3] * state[2] - f[4] * state[3];
        }
    }
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        Opening& opening = openings_[index];
        opening.history = opening.b1 * opening.previous_into - opening.a1 * opening.previous_out;
    }
}

double SurfaceNetwork::glottis_free_pressure() const noexcept {
    double pressure = 2.0 * weighted_[glottis_] / admittance_sum_[glottis_];
    return pressure;
}

double SurfaceNetwork::glottis_impedance() const noexcept {
    double impedance = 1.0 / admittance_sum_[glottis_];
    return impedance;
}

double SurfaceNetwork::scatter(double glottal_flow, double glottal_turbulence, double jet_turbulence) noexcept {
    const std::size_t count = pieces_.size();
    const std::size_t slot = write_index_ & ring_mask_;
    for (std::size_t boundary = 0; boundary < boundary_flow_.size(); ++boundary) {
        boundary_flow_[boundary] = 0.0;
        boundary_admittance_[boundary] = 0.0;
    }
    double wall_flow = 0.0;
    for (std::size_t piece = 0; piece < count; ++piece) {
        // Junction with an instantaneous port of admittance Y whose return is
        // b = (b0 p + h)/(1 + b0) (b0: the filter's direct term, h: its
        // history): p = (W + 2 Y h/(1 + b0)) / (sum Y - 2 Y b0/(1 + b0)).
        double drive = 2.0 * weighted_[piece];
        double total = admittance_sum_[piece];
        std::size_t boundary = piece_boundary_[piece];
        double resistance = 0.0;
        double through = 1.0;
        double turbulence = 0.0;
        if (piece == glottis_) {
            drive += glottal_flow;
        }
        drive += injected_[piece];
        injected_[piece] = 0.0;
        if (boundary != no_parent) {
            // Series resistance and turbulent pressure q on the paths across
            // the boundary: flow in = Y (2a - p - q) / (1 + R Y).
            resistance = boundary_resistance_[boundary];
            through = 1.0 / (1.0 + resistance);
            if (boundary == turbulence_boundary_) {
                turbulence = turbulence_weight_ * jet_turbulence;
            } else if (boundary == fading_boundary_) {
                turbulence = (1.0 - turbulence_weight_) * jet_turbulence;
            }
            total -= resisted_sum_[piece] * (1.0 - through);
            drive -= resisted_sum_[piece] * through * turbulence;
        }
        double port = 0.0;
        double direct = 0.0;
        double history = 0.0;
        std::size_t opening_index = opening_of_piece_[piece];
        double pull = 0.0;
        if (wall_admittance_[piece] > 0.0) {
            // The edges pull with pressure e on the tissue, which then feels
            // p + e: its own waves are the air's plus e/2 each, so the return
            // gains e (b0 - 1)/2 (b0: the filter's direct term).
            port = wall_admittance_[piece];
            direct = wall_filter_[piece][0];
            pull = edge_pressure_[piece];
            history = wall_history_[piece] + 0.5 * pull * (direct - 1.0);
        } else if (opening_index != no_parent) {
            port = openings_[opening_index].admittance;
            direct = openings_[opening_index].b0;
            history = openings_[opening_index].history;
        }
        double pressure = drive / total;
        if (port > 0.0) {
            double share = 1.0 / (1.0 + direct);
            pressure = (drive + 2.0 * port * history * share) / (total - 2.0 * port * direct * share);
            double returned = (direct * pressure + history) * share;
            double into = pressure - returned;
            if (opening_index != no_parent) {
                Opening& opening = openings_[opening_index];
                opening.previous_into = into;
                opening.previous_out = returned;
                opening.flow = port * (into - returned);
            } else {
                wall_flow += port * (into - returned);
                std::array<double, 4>& state = wall_state_[piece];
                state[1] = state[0];
                state[0] = into + 0.5 * pull;
                state[3] = state[2];
                state[2] = returned + 0.5 * pull;
                displacement_[piece] += (into - returned) / (air_.density * air_.sound_speed * sample_rate_);
            }
        }
        pressure_[piece] = pressure;
        // Emission along every path: pressure minus what that path brought in.
        // Turbulence of the glottal jet enters the paths leaving the glottis.
        double source = 0.0;
        if (piece == glottis_) {
            source = glottal_turbulence;
        }
        for (std::size_t position = outgoing_start_[piece]; position < outgoing_start_[piece + 1]; ++position) {
            std::size_t link = outgoing_links_[position];
            std::size_t reverse = link ^ 1u;
            double arriving = arriving_[reverse];
            double emitted = pressure - arriving + source;
            if (link >= first_cross_) {
                emitted = root_admittance_[link - first_cross_] * (pressure + source) - arriving;
            }
            if (link_boundary_[reverse] != no_parent) {
                double admittance = links_[reverse].admittance;
                emitted = (pressure + turbulence - arriving * (1.0 - resistance)) * through;
                boundary_flow_[boundary] += admittance * through * (2.0 * arriving - pressure - turbulence);
                boundary_admittance_[boundary] += admittance;
            }
            history_[links_[link].buffer_offset + slot] = emitted;
        }
    }
    write_index_ = (write_index_ + 1) & ring_mask_;
    wall_flow_ = wall_flow;
    if (turbulence_weight_ < 1.0) {
        turbulence_weight_ = std::min(1.0, turbulence_weight_ + turbulence_fade_);
        if (turbulence_weight_ >= 1.0) {
            fading_boundary_ = no_parent;
        }
    }
    // Each boundary's resistance for the next sample from this sample's flow:
    // R Y = rho |U| K / 2 times the boundary's path admittance, with
    // K = (1/A_narrow - 1/A_wide)^2 when the flow widens across it, else 0.
    for (std::size_t index = 0; index < boundary_flow_.size(); ++index) {
        boundary_smoothed_[index] += smoothing_ * (boundary_flow_[index] - boundary_smoothed_[index]);
        double flow = boundary_smoothed_[index];
        std::size_t before = boundary_upstream_[index];
        std::size_t after = boundary_downstream_[index];
        if (before == no_parent) {
            boundary_resistance_[index] = 0.0;
            continue;
        }
        double inverse_up = 1.0 / geometry_area_[before];
        double inverse_down = after != no_parent ? 1.0 / geometry_area_[after] : 0.0;
        double widening = 0.0;
        if (flow > 0.0 && inverse_down < inverse_up) {
            widening = inverse_up - inverse_down;
        } else if (flow < 0.0 && inverse_up < inverse_down) {
            widening = inverse_down - inverse_up;
        }
        // No speech jet is faster than about 60 m/s (2 kPa of drive). The
        // smoothed flow lags a closing section, so its speed is bounded
        // there; otherwise a sealing lip would see a huge resistance.
        double narrow = 1.0 / std::max(inverse_up, inverse_down);
        double magnitude = std::min(std::fabs(flow), fastest_jet * narrow);
        boundary_resistance_[index] = 0.5 * air_.density * magnitude * widening * widening
            * boundary_admittance_[index];
    }
    double heard = 0.0;
    for (std::size_t index = 0; index < openings_.size(); ++index) {
        Opening& opening = openings_[index];
        heard += (opening.flow - opening.previous_flow) * sample_rate_;
        opening.previous_flow = opening.flow;
    }
    return heard;
}

double SurfaceNetwork::opening_area(std::size_t opening) const noexcept {
    return openings_[opening].area;
}

double SurfaceNetwork::opening_resistance(std::size_t opening) const noexcept {
    return boundary_resistance_[stations_.size() + opening];
}

double SurfaceNetwork::entry_resistance(std::size_t station) const noexcept {
    return boundary_resistance_[station];
}

double SurfaceNetwork::wall_flow() const noexcept {
    return wall_flow_;
}

double SurfaceNetwork::station_pressure(std::size_t station) const noexcept {
    const Station& ring = stations_[station];
    double sum = 0.0;
    for (std::size_t piece = 0; piece < ring.pieces; ++piece) {
        sum += pressure_[ring.first_piece + piece];
    }
    double mean = sum / static_cast<double>(ring.pieces);
    return mean;
}

double SurfaceNetwork::radiated_flow(std::size_t opening) const noexcept {
    return openings_[opening].flow;
}

std::size_t SurfaceNetwork::opening_count() const noexcept {
    return openings_.size();
}

std::size_t SurfaceNetwork::station_count() const noexcept {
    return stations_.size();
}

std::size_t SurfaceNetwork::tube_station(std::size_t tube, std::size_t station) const noexcept {
    std::size_t flat = tube_first_station_[tube] + station;
    return flat;
}

std::size_t SurfaceNetwork::piece_count() const noexcept {
    return pieces_.size();
}

std::size_t SurfaceNetwork::link_count() const noexcept {
    return links_.size();
}

}

namespace ducklips {

void SurfaceNetwork::set_openings(std::span<const double> opening_area) noexcept {
    for (std::size_t flat = 0; flat < stations_.size() && flat < opening_area.size(); ++flat) {
        opening_target_[flat] = beside_[flat] != no_parent ? std::max(opening_area[flat], 0.0) : 0.0;
    }
}

void SurfaceNetwork::add_flow(std::size_t piece, double flow) noexcept {
    injected_[piece] += flow;
}

double SurfaceNetwork::piece_pressure(std::size_t piece) const noexcept {
    return pressure_[piece];
}

std::size_t SurfaceNetwork::first_piece(std::size_t station) const noexcept {
    return stations_[station].first_piece;
}

SurfaceLinearSnapshot SurfaceNetwork::linear_snapshot() const {
    if (pieces_.empty() || glide_left_ != 0) {
        throw std::logic_error("linear snapshot requires prepared, stationary geometry");
    }
    for (std::size_t station = 0; station < stations_.size(); ++station) {
        if (geometry_area_[station] != station_area_now_[station]
            || geometry_first_[station] != station_area_now_[station]
            || opening_area_[station] != opening_now_[station]
            || opening_first_[station] != opening_now_[station]) {
            throw std::logic_error("linear snapshot requires settled geometry; reshape with glide=1 for a frozen probe");
        }
        if (stations_[station].edge_tension != 0.0) {
            throw std::logic_error("linear snapshot does not yet include membrane edge tension");
        }
    }
    for (double resistance : boundary_resistance_) {
        if (resistance != 0.0) {
            throw std::logic_error("linear snapshot requires zero-mean-flow kinetic resistance; reset before probing");
        }
    }
    SurfaceLinearSnapshot snapshot{};
    snapshot.sample_rate = sample_rate_;
    snapshot.pieces = pieces_.size();
    snapshot.glottis = glottis_;
    snapshot.allocated_delay_samples = history_.size();
    snapshot.piece_tubes.reserve(pieces_.size());
    for (const SurfacePiece& piece : pieces_) {
        snapshot.piece_tubes.push_back(stations_[piece.station].tube);
    }
    snapshot.paths.reserve(links_.size());
    for (std::size_t index = 0; index < links_.size(); ++index) {
        const AirLink& link = links_[index];
        SurfaceLinearPath path{};
        path.from = link.from;
        path.to = link.to;
        path.admittance = link.admittance;
        if (index < first_cross_) {
            path.delay = link.whole_delay + 1;
            const std::array<double, 3>& loss = loss_filter_[index];
            path.numerator = {link.allpass * loss[0], loss[0] + link.allpass * loss[1], loss[1]};
            path.denominator = {1.0, link.allpass + loss[2], link.allpass * loss[2]};
        } else {
            double delay = cross_delay_[index - first_cross_];
            path.delay = static_cast<std::size_t>(std::floor(delay));
            double fraction = delay - static_cast<double>(path.delay);
            path.numerator = {1.0 - fraction, fraction, 0.0};
        }
        snapshot.paths.push_back(path);
    }
    for (std::size_t piece = 0; piece < pieces_.size(); ++piece) {
        if (wall_admittance_[piece] > 0.0) {
            const std::array<double, 5>& filter = wall_filter_[piece];
            SurfaceLinearPort port{};
            port.piece = piece;
            port.admittance = wall_admittance_[piece];
            port.numerator = {filter[0], filter[1], filter[2]};
            port.denominator = {1.0, filter[3], filter[4]};
            snapshot.ports.push_back(port);
        }
    }
    for (const Opening& opening : openings_) {
        SurfaceLinearPort port{};
        port.piece = opening.piece;
        port.admittance = opening.admittance;
        port.numerator = {opening.b0, opening.b1, 0.0};
        port.denominator = {1.0, opening.a1, 0.0};
        port.radiation = true;
        snapshot.ports.push_back(port);
    }
    return snapshot;
}

}
