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

# Gallery
<img width="1920" height="1080" alt="Denoise" src="https://github.com/user-attachments/assets/d5f8beb5-69f9-445a-9e25-f9cd8929a5a0" />
Spectral Path Tracing with convergence

<img width="3840" height="2160" alt="4K" src="https://github.com/user-attachments/assets/7196fd2c-9ca6-4d4c-ad76-07f289df3f2f" />
Before including Hero's Spectral Sampling

<img width="1280" height="720" alt="20 00 thinking 1024 spp" src="https://github.com/user-attachments/assets/a67b81c0-8f84-4ac8-815b-b228f8d73242" />
Thinking about Teapots

<img width="1920" height="1080" alt="64 spp 1016" src="https://github.com/user-attachments/assets/f2df6d98-753a-4528-9bc1-2927aeeb6d09" />

<img width="1280" height="720" alt="22 32 intensity" src="https://github.com/user-attachments/assets/e1f648d1-ce2e-4ae1-8960-2666b180307f" />
Photon Mapping with Final Gathering

