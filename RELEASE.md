# Release 0.6.0

## Major Features And Improvements
* Concat operation now supports special modes to concatenate large number of files in parallel when the downstream operation is passed to sort (reported by @ekg). This mode is ~30-fold faster / thread and is parallelizable.
* Introducing the experimental `aggregation` subroutine. Modelled based on the Datashader pipeline, Tomahawk can now natively construct aggregated datasets of billions of points into datasets of pixels that can be plotted using raster operations efficiently.
* Slicing operator for `two` files has been re-instated. Use `-O` for output type (`two` or `ld`) when invoking the `view` subroutine

## Bug Fixes and Other Changes
* Bug fixes
  * Problem when importing variants within 1-base of each other (reported by @ekg)
  * Memory requirement problems when parallel sort merging very large files (reported by @ekg)
  * Removed default missingness filter (from 0.2 to 1.0) to support importing `bcf` with few samples
  * Fixed problem with linked intervals during slicing. Now correctly selects targetted associations based on both mates (reported by @ekg)

# Release 0.5.3

## Minor Features And Improvements
* Default build now uses SSE4.2 only
* Will correctly build `dylib` on OSX
* Added filtering for maxmimum cell counts when viewing `two` files (flags `-1`,`-2`,`-3`,`-4`)

## Bug Fixes and Other Changes
* Bug fixes
   * Viewing a `two` file with `-l` or `-u` now correctly filters output when passed as the sole parameter

# Release 0.5.2

## Minor Features And Improvements
* Now correctly filters `two` files

# Release 0.5.1

## Minor Features And Improvements
* Meta index for faster queries of `two` files now builds correctly
* Interval slicing for sorted `two` files now outputs filtering hits for the from positions only
* Haplotype counts in binary `two` files are now stored as `double` instead of `float` 

# Release 0.5.0

## Breaking Changes
* All changes are breaking because of structural rechanges. No files generated prior to this release are supported.

## Minor Features And Improvements
* Build process now smoother
* Shared C++ library is now built by default

# Release 0.4.0

## Breaking Changes
* All changes are breaking because of structural rechanges. No files generated prior to this release are supported.

## Major Features And Improvements
* Added large number of new FLAG fields to `two` files
* Reintroduced speed mode for calculations when contingency tables are not required

# Release 0.3.3

## Major Features And Improvements
* None

## Bug Fixes and Other Changes
* Bug fixes
   * Viewing sorted two files fixed for edge cases

# Release 0.3.2

## Breaking Changes
* Using old sorted files does not break functionality but have to be resorted

## Major Features And Improvements
* Sorting `two` files not correctly produces two indices used for fast queries
* Viewing sorted two files is significantly faster

## Bug Fixes and Other Changes
* Bug fixes
   * Sort merge (`tomahawk sort -M`) now produces the correct index

# Release 0.3.1

## Major Features And Improvements
* None

## Bug Fixes and Other Changes
* Bug fixes
   * Sort merge (`tomahawk sort -M`) now produces the correct output
   * `view` ABI command now correctly triggers `-h`/`-H` flag

# Release 0.3.0

## Breaking Changes
* All changes are breaking. There is no backwards compatibility from this release!

## Major Features And Improvements
* All major classes are now implemented as STL-like container and are as decoupled as possible
* External sort file handles now buffer a small amount of data to reduce random access lookups
* `two`/`twk` entries are no longer forced to be accessed by unaligned memory addresses
* Both `two` and `twk` are now self-indexing. All external indices are now invalid

## Bug Fixes and Other Changes
* Bug fixes
  * Too many to list here
* Examples updates:
  * Added R script to demonstrate simple plotting using `base`
  * Added several figures demonstrating some keys concepts of Tomahawk