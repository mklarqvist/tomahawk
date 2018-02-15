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