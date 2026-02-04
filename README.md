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
The following features were independently designed, implemented, and were not part of the course curriculum:

### Spectral Rendering
- Full spectral rendering pipeline
- Wavelength sampling strategies
- Spectral → RGB conversion (CIE & Bruton's approximation)

## Iridescence
- Implemented Airy's Reflectance based on the paper titled "A Practical Extension to Microfacet Theory for the
Modeling of Varying Iridescence"

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

# Gallery
<img width="3840" height="2160" alt="09 43 49 4K 8196 SPP MS" src="https://github.com/user-attachments/assets/1b4ab156-7e8a-4ecd-a701-f3e5254732e0" />
Winning Entry for the Teapot Rendering Competition (Phase 1)
Features the Spectral Renderer<br>
<br />

<img width="1280" height="720" alt="sub 2 1024 SPP DENOISE" src="https://github.com/user-attachments/assets/93917cb8-9e07-4200-b0b9-abb05e388975" />
Juror's Choice Award for the Teapot Rendering Competition (Final Phase)
Features Airy's Reflectance on the hexagonal columns and the correct CIE Wavelength to RGB conversion for my Spectral Renderer on the tea pot<br>
<br />
<img width="1280" height="720" alt="20 00 thinking 1024 spp" src="https://github.com/user-attachments/assets/a67b81c0-8f84-4ac8-815b-b228f8d73242" />
Thinking about Teapots<br>
<br />

<img width="1280" height="720" alt="22 32 intensity" src="https://github.com/user-attachments/assets/e1f648d1-ce2e-4ae1-8960-2666b180307f" />
Photon Mapping with Final Gathering

## References
Laurent Belcour, Pascal Barla. A Practical Extension to Microfacet Theory for the Modeling of
Varying Iridescence. ACM Transactions on Graphics, 2017, 36 (4), pp.65. ￿10.1145/3072959.3073620￿.
￿hal-01518344v2￿
