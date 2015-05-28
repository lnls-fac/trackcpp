
#include "interface.h"

Status::type track_elementpass_wrapper (
         const Element& el,
         Pos<double> &orig_pos,
         const Accelerator& accelerator) {
         return track_elementpass (el,
                                   orig_pos,
                                   accelerator);
}

Status::type track_linepass_wrapper(
        const Accelerator &accelerator,
        Pos<double> &orig_pos,
        std::vector< Pos<double> >& pos,
        LinePassArgs& args) {
    return track_linepass(accelerator,
                          orig_pos,
                          pos,
                          args.element_offset,
                          args.lost_plane,
                          args.trajectory);
}

Status::type track_ringpass_wrapper (
        const Accelerator& accelerator,
        Pos<double> &orig_pos,
        std::vector< Pos<double> >& pos,
        RingPassArgs& args) {
    return track_ringpass(accelerator,
                          orig_pos,
                          pos,
                          args.nr_turns,
                          args.lost_turn,
                          args.element_offset,
                          args.lost_plane,
                          args.trajectory);
}
