#pragma once

#include <Message_ProgressRange.hxx>

// A wall-clock budget for one STEP export.
//
// A budget cannot be enforced by checking the clock between operations: a
// single boolean or a single sewing pass is what runs for minutes, and by
// the time it returns the budget has long gone. OCCT's long algorithms
// poll their progress indicator instead, so the budget is expressed as an
// indicator whose UserBreak() turns true once the deadline passes. The
// algorithm then abandons its work and reports failure, and the caller
// takes whichever fallback it would have taken for any other failure.
//
// The budget covers B-rep work only. The mesh fallback it degrades to is
// never cut short, since it is the thing that still produces a file.
//
// It is off by default, and deliberately so: the same model would
// otherwise export exact surfaces on a fast machine and a mesh on a slow
// one, which is a poor property for a file format people diff and archive.
namespace OcctProgress {

// Starts a budget of `seconds` from now; 0 or less clears it. Call once
// per export, and clear() when it ends.
void setBudget(double seconds);
void clear();

// A range to hand to an OCCT algorithm. Returns an unbounded range when no
// budget is set, which is what every algorithm gets by default anyway.
Message_ProgressRange range();

// Did a budget expire during this export? Callers use it to explain a
// failure that is really a deadline.
bool stopped();

}  // namespace OcctProgress
