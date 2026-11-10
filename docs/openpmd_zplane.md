# Position (z = const) snapshots in openPMD

Proposal for the openPMD standard 2.0 (pull request against the `upcoming-2.0.0` branch of
github.com/openPMD/openPMD-standard), by the developers of GEM-PIC-RZ. It follows issue #161 ("Make `time`, `dt`
and `timeUnitSI` optional in basePath", D. Sagan: these do not make sense "when particle positions are recorded at
constant longitudinal position") and the per-particle `time` record of the BeamPhysics extension. This file
explains the proposal and what GEM-PIC-RZ writes now.

## Summary

An iteration (a *snapshot* in 2.0) may hold data recorded at a **fixed position** instead of a fixed time: the
history of everything that crossed a plane, e.g. z = s.
- **Snapshot kind:** the new optional iteration attribute `snapshotType` is `"time"` (default, today's meaning) or
  `"position"`.
- **Plane:** `snapshotAxis` (the axis label normal to the plane, e.g. `"z"`), `snapshotPosition` and
  `snapshotPositionUnitSI` (its position, as `time` and `timeUnitSI`).
- **Time along the axis:** meshes keep their spatial axis along `snapshotAxis`, but its coordinate x is a label
  for the time at which the data crossed the plane,
  `t = time + (snapshotPosition - x) / snapshotAxisVelocity`,
  with `snapshotAxisVelocity` in `snapshotPositionUnitSI / timeUnitSI`. For codes in co-moving variables
  ζ = z − v t this is exactly their native picture: x = ζ + s, t = (s − ζ)/v.
- **Particles:** every particle carries the time at which it crossed the plane in the record `time` (as defined
  by the BeamPhysics extension: absolute time = `time` + `timeOffset`, both per-particle records). Along
  `snapshotAxis` the particle position is the same label as for the meshes, so fields and particles overlay.
- **Compatibility:** all old attributes keep valid values (`time` = the reference time, `dt`, the axis label
  `"z"`), so readers that ignore the new attributes show the data as a t = const snapshot of the co-moving window,
  which is what users of such codes plot anyway.

## Motivation

- **Codes that advance in a longitudinal coordinate, not in time.** GEM-PIC-RZ (Galilean GT2 scheme) advances
  the fields of the window from level to level, a level being the plane z = s; every slice ζ = z − ct of the
  window is known at that plane, i.e. at the time t = (s − ζ)/c. One output step is therefore the time history
  of the plane z = s, not a snapshot at t = s/c. The same holds for FEL and beam-transport codes in the
  coordinates (z, t − z/c), for particle trackers that integrate in s, and for codes in co-moving variables that
  output the window at a fixed longitudinal position of its reference point.
- **Screen diagnostics.** Particle screens (records of the particles that cross a plane, with their crossing
  times), field probes at a plane, streak-camera-like diagnostics of time-domain PIC codes and measured data from
  screens have the same structure.
- **Today** openPMD requires `time`, `dt` and `timeUnitSI` and gives every record a `timeOffset` relative to
  `time`: there is no way to say that the data are at a fixed position and that the time varies along an axis.
  Writers either state a wrong t = const snapshot (GEM-PIC-RZ up to version 1.x of its openPMD output; reported by
  a user comparing with FBPIC) or invent private conventions.
- **Consequences of the wrong claim:** a reader that compares a z = const output with a t = const output of
  another code at "the same time" compares different physics as soon as the pulse evolves while it crosses the
  window (laser depletion, dephasing of a bunch, flying focus): the head of the window is at the plane earlier
  than the tail.

## Specification (text of the pull request)

**New optional attributes of the iteration (`basePath`)**
- `snapshotType` *(string)*: `"time"` (default if absent: today's meaning) or `"position"`.
- If `"position"`, required:
  - `snapshotAxis` *(string)*: the label of the axis normal to the plane (as in `axisLabels` of the meshes and the
    components of `position`), e.g. `"z"`.
  - `snapshotPosition` *(floatX)*: position of the plane along `snapshotAxis`.
  - `snapshotPositionUnitSI` *(float64)*: factor to SI (meters) for `snapshotPosition`.
  - `snapshotAxisVelocity` *(floatX)*: the velocity v that maps the coordinate along `snapshotAxis` to time, in
    `snapshotPositionUnitSI / timeUnitSI`. The value `0` (particles only) means that the coordinate along
    `snapshotAxis` is the true particle position and the time is given only by the `time` record.
- Meaning of the existing attributes in a position snapshot:
  - `time`: the reference time, at which the point x = `snapshotPosition` of the axis crosses the plane.
  - `dt`: as now, the step that led to this snapshot, converted to time (e.g. Δz / v).
  - `timeOffset` of a record: as now, an offset of the record relative to `time` (e.g. a staggered field component
    recorded half a step later).

**Mesh records in a position snapshot**
- The axis `snapshotAxis` stays a regular (or rectilinear) mesh axis with its `gridSpacing`, `gridGlobalOffset`
  and `gridUnitSI`. A data point at the coordinate x (in SI: (gridGlobalOffset + index · gridSpacing + position)
  · gridUnitSI) holds the value at the plane at the time
  `t = time · timeUnitSI + (snapshotPosition · snapshotPositionUnitSI - x) / (snapshotAxisVelocity · snapshotPositionUnitSI / timeUnitSI)`.
- The other axes are spatial coordinates in the plane, as now.

**Particle records in a position snapshot**
- `time` *(scalar record, unitDimension of time)*: time at which the particle crossed the plane, per particle
  (BeamPhysics); absolute time = `time` + `timeOffset` (per-particle record of BeamPhysics, zero if absent).
- `position` + `positionOffset` along `snapshotAxis`: for v ≠ 0 the label x = snapshotPosition − v (t − time),
  as for the meshes; for v = 0 the true position (at the plane: `snapshotPosition` for every particle).
- Particles that are not exactly at the plane (e.g. a code that records the last state before the plane) store
  their true position for v = 0, or the label of their own `time` for v ≠ 0; the `time` record is always the
  time of the stored state.

**Record-level override (optional)**
- The five attributes may also be set on a mesh record or on a particle species group and then take precedence
  over the iteration's: a time-domain PIC code can write t = const meshes and a particle screen in the same
  iteration.

## Reader behaviour

- Readers that do not know the attributes see a t = const snapshot at `time`, with the longitudinal axis
  labelled by the co-moving coordinate: plots and lineouts are correct as functions of the label.
- Readers that know them should
  - label the axis as time (or as the co-moving coordinate ζ = x − snapshotPosition) and show the plane position
    instead of the time in titles;
  - convert the axis to time with the formula above when the user asks for it, and use the `time` record for
    particles;
  - not mix t = const and z = const snapshots in comparisons without warning.
- A **t = const snapshot can be rebuilt** from a series of position snapshots: the slice x of the snapshot at
  plane s holds the state at t(s, x); collecting, for a wanted time t*, the slice from the plane s at which
  t(s, x) = t* (interpolating between neighbouring snapshots) gives the lab-frame picture. This needs output at
  every plane the window crosses within its length L, i.e. about L/Δz snapshots, and is better done in the code
  (see below).

## Example (GEM-PIC-RZ, iteration 4, `h5dump`-style)

```
/data/4                                      # one level: the plane z = s = 160 k0^-1
  time                   = 160.0             # s / c in units 1/omega0: slice zeta = 0 crosses the plane
  dt                     = 40.0
  timeUnitSI             = 4.2441e-16
  snapshotType           = "position"
  snapshotAxis           = "z"
  snapshotPosition       = 160.0             # s in units 1/k0
  snapshotPositionUnitSI = 1.2732e-07
  snapshotAxisVelocity   = 1.0               # c = (1/k0) / (1/omega0)
/data/4/meshes/E                             # thetaMode, axisLabels ["r", "z"]
  gridGlobalOffset       = [0.0, 160.0]      # z label = zeta + s: slice zeta crossed the plane at t = (s - zeta)/c
/data/4/particles/witness/time               float64 [20000], unitSI = 4.2441e-16, (s - zeta)/c per particle
/data/4/particles/witness/position/z         zeta,  positionOffset/z = s: the same label as the meshes
```

## What GEM-PIC-RZ writes now (openPMD 1.1 files)

openPMD 1.1 allows additional attributes and records, so the files stay valid 1.1 files (validator: no errors)
and openpmd-viewer, openPMD-api and VisualPIC read them as before.
- Iteration attributes `snapshotType = "position"`, `snapshotAxis = "z"`, `snapshotPosition` = s,
  `snapshotPositionUnitSI` = 1/k0 in m, `snapshotAxisVelocity = 1` (c); `time` = s/c.
- Meshes: axis label `"z"` = ζ + s (`gridGlobalOffset` = s), as before.
- Particles: the record `time` = (s − ζ)/c (float64), positions as before (`position/z` = ζ,
  `positionOffset/z` = s).
- The root attribute `comment` states the same in words.
- The files do not declare the BeamPhysics extension (it requires the SpeciesType extension and further
  attributes); only its `time` record is used.

Tests: ctest `openpmd` (attributes, particle time against ζ of the native file) and `openpmd_viewer` (all
iterations read with openpmd-viewer, fields and particle moments against the native files).

## Not done: reconstructed t = const output (option)

GEM-PIC-RZ could also write true t = const snapshots: slice ζ is at the lab time t* when the march reaches the
level s = c t* + ζ. Copying each slice at that level (interpolated between the two levels around it) into an
output buffer fills the snapshot after the window has moved by its length L, i.e. after about L/Δz steps;
several snapshots in flight need one buffer each. Such files would be ordinary openPMD files with
`snapshotType = "time"` and the lab coordinate z = ct* + ζ.

## Open questions

- Names: `snapshot…` follows the renaming of iterations to snapshots planned for 2.0 (#148); `iteration…` would
  fit 1.x.
- Planes normal to other axes (x = const screens) are covered by `snapshotAxis`; curved screens are not.
- Whether `time`, `dt`, `timeUnitSI` become optional for position snapshots (#161) or keep the reference-time
  meaning proposed here (better for old readers).
