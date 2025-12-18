# Avohkii-Renderer

A modular path tracer developed for a university rendering course, with
significant independent extensions beyond the course scope.

## Course Scope
The base renderer was implemented as part of a university course and includes:
- RGB path tracing
- Standard BSDFs
- Monte Carlo light transport

(The course framework and provided infrastructure are not included in this repository.)

## Independent Extensions (Beyond Course Scope)
The following features were independently designed and implemented and were
not part of the course curriculum:

### Spectral Rendering
- Full spectral rendering pipeline
- Wavelength sampling strategies
- Spectral → RGB conversion

### Additional Rendering Features
- Environment Map lighting
- Emission and Absorption modeling
- Möller–Trumbore ray–triangle intersection algorithm

These components were implemented for the **Teapot Rendering Competition**.

## Authorship
All functional rendering logic presented in this repository, including all
independent extensions listed above, was implemented by me.

Instructor-provided framework code is not included.

## Build Status
This repository is not standalone and will not build on its own.
