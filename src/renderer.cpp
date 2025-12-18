#include <iostream>

#include <thread>
#include <vector>
#include <atomic>

#include "renderer.h"
#include "scene.h"
#include "objects.h"
#include "materials.h"
#include "lights.h"
#include "rng.h"
#include "photonmap.h"

# define LEGACY_SHADING_API

#define BOUNCES 5
#define PI 3.1415926f

#define BIAS 0.0002f
#define EPSILON 0.00001f

#define SPP_MIN 16
#define SPP_MAX 16
#define DELTA_MAX 0.01f // 1/256 bit size?

Node* firstNode;
PhotonMap apple;
PhotonMap orange;
void NodeTraversal(Node* root, Ray camToPixel, HitInfo& hInfo, int hitSide);

Color WavelengthToRGB(float lambda)// picking a random wavelenght form the visible color spectrum
{
	float r, g, b;


	if (lambda >= 0.38f && lambda < 0.44f)
	{
		r = -(lambda - 0.44f) / (0.44f - 0.38f); //violet to blue
		g = 0.0f;
		b = 1.0f;
	}

	else if (lambda >= 0.44f && lambda < 0.49f) //blue to cyan	
	{
		r = 0.0f;
		g = (lambda - 0.44f) / (0.49f - 0.44f);
		b = 1.0f;
	}

	else if (lambda >= 0.49f && lambda < 0.51f) //cyan to greeeeeen
	{
		r = 0.0f;
		g = 1.0f;
		b = -(lambda - 0.51f) / (0.51f - 0.49f);
	}

	else if (lambda >= 0.51f && lambda < 0.58f) //green to yellow
	{
		r = (lambda - 0.51f) / (0.58f - 0.51f);
		g = 1.0f;
		b = 0.0f;
	}

	else if (lambda >= 0.58f && lambda < 0.64f) //yellow to red
	{
		r = 1.0f;
		g = -(lambda - 0.64f) / (0.64f - 0.58f);
		b = 0.0f;
	}


	else if (lambda >= 0.64f && lambda <= 0.78f) //Just Red
	{
		r = 1.0f;
		g = 0.0f;
		b = 0.0f;
	}

	else //b l a ck can't see shit
	{
		r = g = b = 0.0f;
	}

	float factor;

	if (lambda >= 0.38f && lambda < 0.42f)//shitty vision
		factor = 0.3f + 0.7f * (lambda - 0.38f) / (0.42f - 0.38f);

	else if (lambda >= 0.42f && lambda < 0.70f)  //2020 baby
		factor = 1.0f;

	else if (lambda >= 0.70f && lambda <= 0.78f)//shitty vision
		factor = 0.3f + 0.7f * (0.78f - lambda) / (0.78f - 0.70f);

	else//can't see shit again
		factor = 0.0f;

	return Color(r * factor, g * factor, b * factor);
}


float PowerHeuristic(float pdfA, float pdfB)
{
	float a = pdfA * pdfA;
	float b = pdfB * pdfB;
	return a / (a + b);
}


float BalanceHeuristic(float pdfA, float pdfB)
{
	float a = pdfA;
	float b = pdfB;
	return a / (a + b);
}

bool ShadowTraverse(Node* root, const Ray& ray, float t_max)
{
	bool block = false;

	for (int h = 0; h < root->GetNumChild(); h++)
	{
		Node* child = root->GetChild(h);
		Ray transRay = child->ToNodeCoords(ray);
		bool block = ShadowTraverse(child, transRay, t_max);
		HitInfo hInfo;

		if (child->GetNodeObj() != nullptr)
		{
			if (child->GetNodeObj()->IntersectRay(transRay, hInfo, 1))
			{
				if (hInfo.z < t_max && hInfo.z > BIAS)
				{
					return true;
				}

				else
				{
					continue;
				}
			}
		}
	}

	return block;
}

bool ShadeInfoExtend::CanBounce() const
{
	if (CurrentBounce() < BOUNCES)
	{
		return true;
	}

	return false;
}


float ShadeInfoExtend::TraceShadowRay(Ray const& ray, float t_max) const
{
	if (ShadowTraverse(firstNode, ray, t_max))
	{
		return 0.0;
	}

	else
	{
		return 1.0;
	}
}

Color ShadeInfoExtend::TraceSecondaryRay(Ray const& ray, float& dist, bool reflection) const
{
	if (CanBounce())
	{
		HitInfo reHitInfo;
		reHitInfo.node = firstNode;
		reHitInfo.z = std::numeric_limits<float>::max();
		NodeTraversal(firstNode, ray, reHitInfo, 0);
		ShadeInfoExtend reInfo = *this;
		reInfo.IncrementBounce();

		reInfo.SetHit(ray, reHitInfo);

		if (reHitInfo.node->GetNodeObj() != nullptr) //maybe check closest.front so we don't reflect back hits maybe?
		{
			return reHitInfo.node->GetMaterial()->Shade(reInfo);
		}
	}

	return EvalEnvironment(ray.dir); //it should return environment color for textures
}

Vec3f CosWeightedHemisphereSampling(Vec3f N, ShadeInfo const& pInfo)
{
	float xHemi = pInfo.RandomFloat();
	float phiHemi = pInfo.RandomFloat() * 2.0f * PI;

	float hemiCosTheta = Sqrt(1.0f - xHemi);// 1.0f - xHemi;
	float hemiSinTheta = Sqrt(xHemi);// sqrtf(1.0f - hemiCosTheta * hemiCosTheta);

	Vec3f x, y;
	N.GetOrthonormals(x, y);
	Vec3f hemiSample = N * hemiCosTheta + x * hemiSinTheta * cosf(phiHemi) + y * hemiSinTheta * sinf(phiHemi);
	hemiSample.Normalize();

	return hemiSample;
}

void PhotonBouncing(Node* root, Ray photonRay, HitInfo& hInfo)
{
	if (root->GetNumChild() <= 0)
	{
		return;
	}

	for (int j = 0; j < root->GetNumChild(); j++)
	{
		Node* child = root->GetChild(j);
		Ray rLocal = child->ToNodeCoords(photonRay);
		PhotonBouncing(child, rLocal, hInfo);

		if (child->GetNodeObj() != nullptr)
		{
			if (child->GetNodeObj()->IntersectRay(rLocal, hInfo, 0))
			{
				hInfo.node = child;
				hInfo.node->FromNodeCoords(hInfo);
			}
		}

	}

	for (int h = 0; h < root->GetNumChild(); h++)
	{
		if (root->GetChild(h) == hInfo.node)
		{
			root->FromNodeCoords(hInfo);
			break;
		}
	}

	return;
}

void PopulatePhotonMap(Ray photonRay, Color power, HitInfo& hInfo, ShadeInfoExtend& pInfo)
{
	apple.AddPhoton(hInfo.p, -photonRay.dir, power);
}

void PopulateCausticsMap(Ray photonRay, Color power, HitInfo& hInfo, ShadeInfoExtend& pInfo)
{
	orange.AddPhoton(hInfo.p, -photonRay.dir, power);
}

PhotonMap const* RayTracer::GetPhotonMap() const
{
	return &apple;
}


PhotonMap const* RayTracer::GetCausticsMap() const
{
	return &orange;
}

void RayTracer::BeginRender()
{
	firstNode = &scene.rootNode;
	float aspectRatio = (float)camera.imgWidth / camera.imgHeight;

	float l = 1.0f;// camera.focaldist; //distance from camera to image display.

	float h = 2.0f * l * tan((camera.fov / 2.0f) * (PI / 180));
	float w = h * aspectRatio;
	int width = this->renderImage.GetWidth(); //IMAGE HEIGHT
	int height = this->renderImage.GetHeight(); //IMAGE WIDTH
	float s = w / (float)camera.imgWidth; //pixel dimension

	int camHeight = camera.imgHeight;
	int camWidth = camera.imgWidth;
	float camH = 2.0 * l * tan(camera.fov / 2.0f) * (PI / 180);
	float camW = camH * (camWidth / camHeight);

	Color24* pixel = this->renderImage.GetPixels();
	float* zBuffer = this->renderImage.GetZBuffer();
	int* sampleCount = this->renderImage.GetSampleCount();

	//=============== Generating Camera Rays in WWWOOOORLD SPAAACE ===============//
	Vec3f cameraZ = -camera.dir;
	Vec3f cameraY = camera.up;
	Vec3f cameraX = cameraY.Cross(cameraZ); //vector orthogonal to Z and Y
	Vec3f topLeft = camera.pos - (w / 2.0) * cameraX + (h / 2.0) * cameraY - l * cameraZ; //literally the top left of ur image


	//=============== OVER-CONFIDENCE ===============//
	float CI80[128] = {
	3.078, 1.886, 1.638, 1.533, 1.476, 1.440, 1.415, 1.397, 1.383, 1.372,
	1.363, 1.356, 1.350, 1.345, 1.341, 1.337, 1.333, 1.330, 1.328, 1.325,
	1.323, 1.321, 1.319, 1.318, 1.316, 1.315, 1.314, 1.313, 1.311, 1.310,
	1.308, 1.307, 1.305, 1.304, 1.303, 1.302, 1.301, 1.300, 1.299, 1.298,
	1.297, 1.296, 1.295, 1.294, 1.293, 1.292, 1.291, 1.290, 1.289, 1.288,
	1.287, 1.286, 1.285, 1.284, 1.283, 1.282, 1.281, 1.280, 1.279, 1.278,
	1.277, 1.276, 1.275, 1.274, 1.273, 1.272, 1.271, 1.270, 1.269, 1.268,
	1.267, 1.266, 1.265, 1.264, 1.263, 1.262, 1.261, 1.260, 1.259, 1.258,
	1.257, 1.256, 1.255, 1.254, 1.253, 1.252, 1.251, 1.250, 1.249, 1.248,
	1.247, 1.246, 1.245, 1.244, 1.243, 1.242, 1.241, 1.240, 1.239, 1.238,
	1.237, 1.236, 1.235, 1.234, 1.233, 1.232, 1.231, 1.230, 1.229, 1.228,
	1.227, 1.226, 1.225, 1.224, 1.223, 1.222, 1.221, 1.220, 1.219, 1.218,
	1.217, 1.216, 1.215, 1.214, 1.213, 1.212, 1.211, 1.210 };

	float CI95[128] = {
	12.706f, 4.303f, 3.182f, 2.776f, 2.571f, 2.447f, 2.365f, 2.306f, 2.262f, 2.228f,
	2.201f, 2.179f, 2.160f, 2.145f, 2.131f, 2.120f, 2.110f, 2.101f, 2.093f, 2.086f,
	2.080f, 2.074f, 2.069f, 2.064f, 2.060f, 2.056f, 2.052f, 2.048f, 2.045f, 2.042f,
	2.039f, 2.036f, 2.033f, 2.031f, 2.029f, 2.027f, 2.025f, 2.023f, 2.021f, 2.020f,
	2.018f, 2.017f, 2.016f, 2.015f, 2.014f, 2.013f, 2.012f, 2.011f, 2.010f, 2.009f,
	2.008f, 2.007f, 2.006f, 2.005f, 2.004f, 2.003f, 2.002f, 2.001f, 2.000f, 1.999f,
	1.998f, 1.997f, 1.996f, 1.995f, 1.994f, 1.993f, 1.992f, 1.991f, 1.990f, 1.989f,
	1.988f, 1.987f, 1.986f, 1.985f, 1.984f, 1.983f, 1.982f, 1.981f, 1.980f, 1.979f,
	1.978f, 1.977f, 1.976f, 1.975f, 1.974f, 1.973f, 1.972f, 1.971f, 1.970f, 1.969f,
	1.968f, 1.967f, 1.966f, 1.965f, 1.964f, 1.963f, 1.962f, 1.961f, 1.960f, 1.959f,
	1.958f, 1.957f, 1.956f, 1.955f, 1.954f, 1.953f, 1.952f, 1.951f, 1.950f, 1.949f,
	1.948f, 1.947f, 1.946f, 1.945f, 1.944f, 1.943f, 1.942f, 1.941f, 1.940f, 1.939f,
	1.938f, 1.937f, 1.936f, 1.935f, 1.934f, 1.933f, 1.932f, 1.931f };

	float CI99[128] = {
	63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
	3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
	2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750,
	2.744, 2.738, 2.733, 2.728, 2.724, 2.719, 2.715, 2.712, 2.708, 2.704,
	2.700, 2.696, 2.692, 2.688, 2.684, 2.680, 2.676, 2.672, 2.668, 2.664,
	2.660, 2.656, 2.652, 2.648, 2.644, 2.640, 2.636, 2.632, 2.628, 2.624,
	2.620, 2.616, 2.612, 2.608, 2.604, 2.600, 2.596, 2.592, 2.588, 2.584,
	2.580, 2.576, 2.572, 2.568, 2.564, 2.560, 2.556, 2.552, 2.548, 2.544,
	2.540, 2.536, 2.532, 2.528, 2.524, 2.520, 2.516, 2.512, 2.508, 2.504,
	2.500, 2.496, 2.492, 2.488, 2.484, 2.480, 2.476, 2.472, 2.468, 2.464,
	2.460, 2.456, 2.452, 2.448, 2.444, 2.440, 2.436, 2.432, 2.428, 2.424,
	2.420, 2.416, 2.412, 2.408, 2.404, 2.400, 2.396, 2.392, 2.388, 2.384,
	2.380, 2.376, 2.372, 2.368, 2.364, 2.360, 2.356, 2.352 };

	//++++++++++++ BUCKET SET UP ++++++++++++//
	//+++++++++ FOR MULTI-THREADING +++++++++//
	const int TILE_SIZE = 32;
	int numTilesX = (width + TILE_SIZE - 1) / TILE_SIZE;
	int numTilesY = (height + TILE_SIZE - 1) / TILE_SIZE;
	int totalTiles = numTilesX * numTilesY;

	std::atomic<int> tileCounter(0);
	std::atomic<int> pixelsRendered(0);

	//========== Get Thread Count ==========//
	unsigned int numThreads = std::thread::hardware_concurrency();
	if (numThreads == 0) numThreads = 4; //=== Fallback ===//
	std::vector<std::thread> threads;
	//++++++++++++++++++++++++++++++++++++++//
	//++++++++++++++++++++++++++++++++++++++//

	//===================== PHOTON MAPPING NERDS =====================//


	RNG random;
	ShadeInfoExtend pInfo(scene.lights, scene.environment, random);
	int maxDepth = 0;
	int photonsEmitted = 1000000;

	apple.Resize(photonsEmitted * 10);
	orange.Resize(photonsEmitted * 2);

	for (int l = 0; l < pInfo.NumLights(); l++)
	{
		if (pInfo.GetLight(l)->IsPhotonSource())
		{
			for (int p = 0; p < photonsEmitted; p++)
			{
				HitInfo pHitInfo;
				pHitInfo.node = firstNode;
				pHitInfo.z = std::numeric_limits<float>::max();

				//=== Uniformly Sampling a Sphere ===/
				Box lightBox = pInfo.GetLight(l)->GetBoundBox();
				float radius = (fabs(lightBox.pmin.y - lightBox.pmax.y)) / 2.0f;
				Vec3f center = (lightBox.pmin + lightBox.pmax) / 2.0f;

				float xSphere = pInfo.RandomFloat();
				float phiSphere = pInfo.RandomFloat() * 2.0f * PI;
				float cosTheta = 1.0f - 2.0f * xSphere;
				float sinTheta = Sqrt(1.0f - (cosTheta * cosTheta));

				Vec3f randLightDir = Vec3f(sinTheta * cosf(phiSphere), sinTheta * sinf(phiSphere), cosTheta);
				randLightDir.Normalize();

				Color power = (pInfo.GetLight(l)->Intensity()) * 4.0f * radius * radius * PI / (photonsEmitted * 10.0f);

				//=== cos weighted Sampling a Hemisphere ===//
				Vec3f hemiSample = CosWeightedHemisphereSampling(randLightDir, pInfo);

				//=== Create Photon Ray (and calculate the ray origin point) ===//
				Ray photonRay;
				photonRay.p = center + randLightDir * radius;
				photonRay.dir = hemiSample;
				photonRay.dir.Normalize();


				//====================SPECTRAL STUFF====================//
				pInfo.SetPathLambda(0.38f + pInfo.RandomFloat() * 0.40f); //WAVELENGTH THINGS makes sure I'm setting them between the visible spectra 0.38 - 0.78
				//====================SPECTRAL STUFF====================//

				PhotonBouncing(firstNode, photonRay, pHitInfo);
				pInfo.SetHit(photonRay, pHitInfo);

				//Alright for pInfo, pHitInfo, and all my samples filled.
				DirSampler::Info si;
				bool fromSpecular = false;

				for (int depth = 0; depth < maxDepth; depth++)
				{
					//doing this to skip the first hit AKA direct lighting.

					if (pHitInfo.node->GetNodeObj() == nullptr || pHitInfo.z >= std::numeric_limits<float>::max())
					{
						break;
					}

					if (!pHitInfo.node->GetMaterial()->GenerateSample(pInfo, photonRay.dir, si))
					{
						break;
					}

					if (depth < 0) //Change this to get direct lighting
					{
						power = power * si.mult / si.prob;
						PopulatePhotonMap(photonRay, power, pHitInfo, pInfo);
					}


					if (si.lobe == DirSampler::DIFFUSE)
					{
						if (depth == 0)
						{
							fromSpecular = false;
						}

						else if (!fromSpecular)
						{
							PopulatePhotonMap(photonRay, power, pHitInfo, pInfo);
						}

						else
						{
							fromSpecular = false;
							PopulateCausticsMap(photonRay, power, pHitInfo, pInfo);
							break;
						}
					}

					else if (si.lobe == DirSampler::SPECULAR)
					{
						fromSpecular = true;
					}

					else if (si.lobe == DirSampler::TRANSMISSION)
					{
						fromSpecular = true;
					}

					else
					{
						break;
					}

					power = power * si.mult / si.prob;
					pHitInfo.N.Normalize();
					photonRay.p = pInfo.P() + pInfo.N() * EPSILON;
					pHitInfo.z = std::numeric_limits<float>::max();
					pHitInfo.node = firstNode;
					PhotonBouncing(firstNode, photonRay, pHitInfo);
					pInfo.SetHit(photonRay, pHitInfo);
				}
			}
		}
	}

	std::cout << apple.NumPhotons() << std::endl;
	std::cout << orange.NumPhotons() << std::endl;
	apple.PrepareForIrradianceEstimation();
	orange.PrepareForIrradianceEstimation();


	auto renderTile = [&]() {
	//================== CRITICAL: Each thread needs its own RNG instance ==================//
		RNG threadRNG;
		ShadeInfoExtend sInfo(scene.lights, scene.environment, threadRNG);
		sInfo.rootNode = firstNode;

		while (true) {
			int tileIdx = tileCounter.fetch_add(1);
			if (tileIdx >= totalTiles) break;

			int tileX = tileIdx % numTilesX;
			int tileY = tileIdx / numTilesX;

			int startX = tileX * TILE_SIZE;
			int startY = tileY * TILE_SIZE;
			int endX = std::min(startX + TILE_SIZE, width);
			int endY = std::min(startY + TILE_SIZE, height);


			for (int i = startX; i < endX; i++)
			{
				for (int j = startY; j < endY; j++)
				{
					Color cum = Color(0.0f, 0.0f, 0.0f);
					Color cumSq = Color(0.0f, 0.0f, 0.0f);

					Color wMean = Color(0.0f, 0.0f, 0.0f);
					Color wVar = Color(0.0f, 0.0f, 0.0f);

					int index = j * width + i;
					float exit = SPP_MAX;

					//================== 4-D Low-Discrepancy Adaptive Sampling for AA and DoF ==================//
					for (int n = 0; n < SPP_MAX; n++)
					{
						//====================SPECTRAL STUFF====================//
						//====================SPECTRAL STUFF====================//
						float pathLambda = 0.38f + Halton(n, 11) * 0.40f;
						sInfo.SetPathLambda(pathLambda);
						//====================SPECTRAL STUFF====================//
						//====================SPECTRAL STUFF====================//

						exit--;
						Color current = Color(0.0f, 0.0f, 0.0f);

						float rngX = Halton(n + (index * SPP_MAX), 2);
						float rngY = Halton(n + (index * SPP_MAX), 3);

						float rngT = Halton(n + (index * SPP_MAX), 5) * 2.0f * PI;
						float rngR = Halton(n + (index * SPP_MAX), 7);
						Vec2f inAperture = Vec2f(cos(rngT), sin(rngT));
						inAperture *= Sqrt(rngR) * camera.dof;

						HitInfo hInfo;
						hInfo.node = firstNode;
						hInfo.front;
						hInfo.z = std::numeric_limits<float>::max();
						float closestZ = std::numeric_limits<float>::max();
						Vec3f pixelPos = topLeft + (i + rngX) * s * cameraX - (j + rngY) * s * cameraY; //This is the position of the CENTER of your pixel

						Ray cameraRay;
						cameraRay.dir = Vec3f(pixelPos - (camera.pos + cameraX * inAperture.x + cameraY * inAperture.y)); // This gets you a vector pointing FROM camera TO the pixel center. Going OUTWARDS
						cameraRay.dir.Normalize();
						cameraRay.p = camera.pos + cameraX * inAperture.x + cameraY * inAperture.y; // They're ALL gonna be set at the camera pos
						NodeTraversal(firstNode, cameraRay, hInfo, 0);

						float u = -(camW / 2.0f) + (camW * (i + 0.5f) / camWidth);
						float v = (camH / 2.0f) - (camH * (j + 0.5f) / camHeight);

						//================== Rendering Renderable Lights Lol ==================//
						for (int m = 0; m < sInfo.NumLights(); m++)
						{
							if (sInfo.GetLight(m)->IsRenderable())
							{
								if (sInfo.GetLight(m)->IntersectRay(cameraRay, hInfo, 0))
								{
									if (hInfo.z < closestZ)
									{
										hInfo.light = true;
										closestZ = hInfo.z;
										sInfo.SetHit(cameraRay, hInfo);
										zBuffer[index] = hInfo.z;
										current = sInfo.GetLight(m)->Radiance(sInfo);
										cum += current;
 									}
								}
							}
						}

						//================== THE ACTUAL SHADE FUNCTION PART ==================//
						//======================== PATHTRACING NOW!!! ========================//
						if (hInfo.node->GetNodeObj() != nullptr)// && !hInfo.light)
						{
							Color throughput = Color(1.0f, 1.0f, 1.0f); //Starts white and gets spread out since it should be energy conserving
							Color L = Color(0.0f, 0.0f, 0.0f); //Light for MIS i think

							DirSampler::Info si; //si for materials
							DirSampler::Info siLight;// si for lights
							DirSampler::Info siPrev;

							Ray wi = cameraRay; //incoming light "which would be the next direction since light is bouncing off it to eye
							Ray wo = cameraRay; //outgoing direction towrads the eye
							Vec3f Li;
							sInfo.SetHit(cameraRay, hInfo);
							sInfo.rootNode = firstNode;
							zBuffer[index] = hInfo.z;
							bool hitLight = false;
							bool firstDiffuse = true;

							for (int p = 0; p < 16; p++)
							{
								float fogT = std::min(wo.dir.Length() * 50.0f, 10.0f);
								throughput *= expf(-0.03f * fogT);

								wo = wi;
								sInfo.SetHit(wo, hInfo);

								if (sInfo.GetNode() == nullptr || sInfo.GetNode()->GetMaterial() == nullptr)
								{
									break;
								}
								
								if (std::strcmp(sInfo.GetNode()->GetMaterial()->GetName(),"env") == 0)
								{
									Color envColor = Color(0.0f, 0.0f, 0.0f);
									int numSamples = 16; // Adjust for more/less blur
									float blurAngle = 0.0f; // Adjust blur radius (radians)

									for (int e = 0; e < numSamples; e++)
									{
										//Offset Random Direction
										float theta = sInfo.RandomFloat() * 2.0f * PI;
										float phi = sInfo.RandomFloat() * blurAngle;

										// Create orthonormal basis around wo.dir
										Vec3f u, v;
										wo.dir.GetOrthonormals(u, v);

										// Offset again in the cone or w/e
										Vec3f offsetDir = wo.dir * cosf(phi) +
											(u * cosf(theta) + v * sinf(theta)) * sinf(phi);
										offsetDir.Normalize();

										envColor += sInfo.EvalEnvironment(offsetDir);
										current = envColor;
									}

									L += throughput * (envColor / numSamples);
									break;

									//L += throughput * sInfo.EvalEnvironment(wo.dir);
									//break;
								}

								if (sInfo.GetNode()->GetMaterial()->emission() != Color(0.0))
								{
									L += throughput * sInfo.GetNode()->GetMaterial()->emission();
							
								}

								for (int l = 0; l < sInfo.NumLights(); l++)
								{
									if (sInfo.GetLight(l)->IntersectRay(wo, hInfo, 0))
									{
										hitLight = true;
										hInfo.light = true;

										if (p == 0)
										{
											L += throughput * sInfo.GetLight(l)->Radiance(sInfo);// *PowerHeuristic(si.prob, LIGHT_IS.prob) / si.prob;
											current = L;
										}

										else if (p > 0)
										{
											DirSampler::Info LIGHT_IS;
											float coso = std::max(0.0f, sInfo.N().Dot(wo.dir));
											sInfo.GetLight(l)->GetSampleInfo(sInfo, wo.dir, LIGHT_IS);
											
											L += throughput * sInfo.GetLight(l)->Radiance(sInfo) * PowerHeuristic(si.prob, LIGHT_IS.prob) / si.prob;
											current = L;
										}
									}

									if (hitLight)
									{
										break;
									}
								}

								if (hitLight)
								{
									break;
								}

								if (sInfo.GetNode()->GetMaterial()->GenerateSample(sInfo, wi.dir, si))
								{
									hInfo.N.Normalize();

									if (si.lobe == DirSampler::DIFFUSE)
									{
										if (firstDiffuse)
										{
											Color irradCaustics;
											Vec3f dirCaustics;
											Color caustics = Color(0.0, 0.0, 0.0);
											orange.EstimateIrradiance<20, PHOTONMAP_FILTER_QUADRATIC>(irradCaustics, dirCaustics, 0.01f, sInfo.P(), sInfo.N(), 0.95f);
											caustics = irradCaustics * throughput * si.mult;// *-dirCaustics.Dot(sInfo.N());
											L += caustics;
											firstDiffuse = false;
										}

										for (int l = 0; l < sInfo.NumLights(); l++)
										{
											if (sInfo.GetLight(l)->GenerateSample(sInfo, Li, siLight))
											{
												float cosLi = std::max(EPSILON, sInfo.N().Dot(Li));
												DirSampler::Info BSDF_IS;
												sInfo.GetNode()->GetMaterial()->GetSampleInfo(sInfo, Li, BSDF_IS);

												L += throughput * si.mult * cosLi * siLight.mult * PowerHeuristic(siLight.prob, BSDF_IS.prob) / siLight.prob;
												current = L;
											}
										}
									}

									throughput *= (si.mult / si.prob);

									if (si.lobe == DirSampler::DIFFUSE)
									{
										throughput *= sInfo.N().Dot(wi.dir);
									}

									wi.p = sInfo.P() + sInfo.N() * EPSILON;
									hInfo.z = std::numeric_limits<float>::max();
									hInfo.node = firstNode;
									siPrev = si;
									NodeTraversal(firstNode, wi, hInfo, 0);
								}		

								else
								{
									break;
								}
								//current = hInfo.node->GetMaterial()->Shade(sInfo);


							}

							//ONLY THE SPECTRAL STUFF WAS ADDED LITERALLY


							Color spectralRGB = WavelengthToRGB(sInfo.GetPathLambda());
							L *= spectralRGB * (1.0f / 0.40f);
							current = L;
							cum += L;
						}

						/*else if (!hInfo.light)
						{
							current = scene.background.Eval(Vec3f(0.5f * (u + 1.0), 0.5f * (v + 1.0), 0.0f));
							//current = scene.background.Eval(Vec3f((i + 0.5f) / camWidth, (j + 0.5f) / camHeight, 0.0f));
							cum += current;
							zBuffer[index] = BIGFLOAT;
						}*/

						//=============== Incrementally Updating the Sample Mean and Variance ===============//
						//===================== This is mostly for Anti-Aliasing but w.e =====================//

						int num = n + 1;
						Color M1 = wMean;
						wMean += (current - wMean) / (num);
						Color M2 = wMean;
						wVar += (current - M1) * (current - M2);

						if (n >= SPP_MIN - 1)
						{
							float t;

							if (n >= 1024)
							{
								t = 1.930f;
							}

							else
							{
								t = CI95[num - 1];
							}

							float SDr = Sqrt(wVar.r / (num - 1));
							float SDg = Sqrt(wVar.g / (num - 1));
							float SDb = Sqrt(wVar.b / (num - 1));

							float apple = t / Sqrt(num);
							float CIr = apple * SDr;
							float CIg = apple * SDg;
							float CIb = apple * SDb;

							sampleCount[index] = num;

							if (CIr <= DELTA_MAX && CIg <= DELTA_MAX && CIb <= DELTA_MAX)
							{
								n = SPP_MAX + 1;
							}
						}
					}

					Color final = cum / (SPP_MAX - exit);

					//=== Built in Gamma Correction ===//

					if (camera.sRGB)
					{
						final = final.Linear2sRGB();
					}

					pixel[index] = Color24(final);
					pixelsRendered.fetch_add(1);
					this->renderImage.IncrementNumRenderPixel(1);
				}
			}
		}
	};

	//=== Launch Threads ===//
	for (unsigned int t = 0; t < numThreads; t++) 
	{
		threads.emplace_back(renderTile);
	}

	//=== Wait for all threads to complete ===
	for (auto& thread : threads)
	{
		thread.join();
	}

	this->renderImage.ComputeZBufferImage();
	this->renderImage.ComputeSampleCountImage();
	this->renderImage.SaveZImage("C:\\Users\\aznbo\\Desktop\\ballz.png");
	this->renderImage.SaveSampleCountImage("C:\\Users\\aznbo\\Desktop\\sampleball.png");
	this->renderImage.SaveImage("C:\\Users\\aznbo\\Desktop\\balls.png");
}

void RayTracer::StopRender()
{

}