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

### Iridescence
- Implemented Airy's Reflectance based on the paper titled "A Practical Extension to Microfacet Theory for the
Modeling of Varying Iridescence"

### Additional Rendering Features
- Environment Map lighting
- Emission and Absorption modeling
- Möller–Trumbore ray–triangle intersection algorithm

These components were implemented for the **Teapot Rendering Competition**.

## Performance
- Optimized render times from hours to seconds (3,200x speedup) by implementing multi-threading and BVH (Bounding Volume Hierarchy) acceleration structures in complex scenes.
- Implemented Low-Discrepancy Monte-Carlo Sampling (Halton/Sobol Sequences) in Anti-Aliasing, Depth of Field, and Soft Shadows significantly reducing noise and render time.

## Authorship
All functional rendering logic presented in this repository, including all
independent extensions listed above, was implemented by me.

Instructor-provided framework code is not included.

## Build Status
This repository is not standalone and will not build on its own.

## Gallery
<table>
  <tr>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/1b4ab156-7e8a-4ecd-a701-f3e5254732e0" width="100%"/>
      <p><b>Winning entry — Teapot Rendering Competition, Phase 1</b><br/>
      <sub>Spectral renderer, Bruton's approximation for wavelength → RGB</sub></p>
    </td>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/93917cb8-9e07-4200-b0b9-abb05e388975" width="100%"/>
      <p><b>Juror's Choice Award — Final Phase</b><br/>
      <sub>Airy's Reflectance on the hexagonal columns, correct CIE conversion</sub></p>
    </td>
  </tr>
  <tr>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/a67b81c0-8f84-4ac8-815b-b228f8d73242" width="100%"/>
      <p><b>Thinking about teapots</b></p>
    </td>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/2b7d073b-0735-4998-aec0-8a67e4203b5a" />
" width="100%"/>
      <p><b>First iteration of competition entry before including Spectral Rendering. Thought it'd be real funny to have the Creation of Adam but with the teapot in the middle. The idea stuck with me.</b></p>
    </td>
  </tr
      <tr>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/97fddff6-fde6-4a08-b3c8-b4769925d73b" />
" width="100%"/>
      <p><b>Testing Airy's Reflectance on white material </b></p>
    </td>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/e1f648d1-ce2e-4ae1-8960-2666b180307f" width="100%"/>
      <p><b>Photon mapping with final gathering</b></p>
    </td>
  </tr>
</table>

## References
Laurent Belcour, Pascal Barla. A Practical Extension to Microfacet Theory for the Modeling of
Varying Iridescence. ACM Transactions on Graphics, 2017, 36 (4), pp.65. 10.1145/3072959.3073620.
hal-01518344v2

Wilkie, A., Nawaz, S., Droske, M., Weidlich, A., & Hanika, J. (2014). Hero Wavelength Spectral Sampling. 
Computer Graphics Forum, 33(4), 123–131. https://doi.org/10.1111/cgf.12419
