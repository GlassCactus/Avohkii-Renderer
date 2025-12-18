#pragma once

#include "materials.h"
#include "renderer.h"
#include "photonmap.h"
#include <iostream>

#define PI 3.1415926f
#define BIAS 0.0002f
#define EPSILON 0.00001f
#define GI_SAMPLES 1
#define GI_MAX 1

extern PhotonMap apple;
extern PhotonMap orange;
extern Node* firstNode; //A global sin
void NodeTraversal(Node* root, Ray camToPixel, HitInfo& hInfo, int hitSide);


float GetDispersiveIOR(float baseIOR, float lambda, float dispersionStrength)
{
	// Cauchy equation is where IOR(lambda) = A + B/lambda^2
	// lambda is in micro meters

	float B = dispersionStrength; // Typically 0.01 - 0.05 for glass but i can go HIGHER
	return baseIOR + B / (lambda * lambda);
}


Vec3f ReflectRay(SamplerInfo const& pInfo, Vec3f N)
{
	Vec3f reflectRay;
	reflectRay = 2.0f * (pInfo.V().Dot(N)) * N - pInfo.V(); //=== CALCULATING THE REFLECTED RAY ===//
	reflectRay.Normalize();

	return reflectRay;
}

Vec3f CosWeightedHemisphereSampling(SamplerInfo const& pInfo, Vec3f N)
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

Color MtlPhong::Shade(ShadeInfo const& sInfo) const //Lol don't need this anymore
{
	return Color(1.0);
}

bool MtlPhong::GenerateSample(SamplerInfo const& sInfo, Vec3f& dir, Info& si) const
{
	return true;
}

void MtlPhong::GetSampleInfo(SamplerInfo const& sInfo, Vec3f const& dir, Info& si) const
{
	return;
}

Color MtlBlinn::Shade(ShadeInfo const& sInfo) const
{
	Color blinnPhong(0.0, 0.0, 0.0);
	Color reflectColor(0.0, 0.0, 0.0);
	Color refractColor(0.0, 0.0, 0.0);
	Color GI(0.0, 0.0, 0.0);
	Vec3f N = sInfo.N();
	N.Normalize();
	Vec3f V = sInfo.V();
	V.Normalize();
	Vec3f P = sInfo.P();

	Color diffColor = sInfo.Eval(diffuse);
	Color specColor = sInfo.Eval(specular);

	Color irrad;
	Vec3f direction;
	Color indirect = Color(0.0, 0.0, 0.0);

	Color irradCaustics;
	Vec3f dirCaustics;
	Color caustics = Color(0.0, 0.0, 0.0);

	//=== GLOSS FLOSS ===//
	// Also sampling a random normal based on the current normal and using it for reflections and refractions
	float gloss = sInfo.Eval(glossiness);
	float xRand = sInfo.RandomFloat();
	float phi = sInfo.RandomFloat() * 2.0f * PI;
	float cosTheta = powf(1.0f - xRand, 1.0f / (gloss + 1.0f));
	float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
	Vec3f u;
	Vec3f v;
	N.GetOrthonormals(u, v);
	Vec3f microfacetN = N * cosTheta + u * sinTheta * cosf(phi) + v * sinTheta * sinf(phi);
	microfacetN.Normalize();
	Vec3f halfDir;

	if (!sInfo.IsFront())
	{
		N = -N;
		microfacetN = -microfacetN;
	}

	for (int i = 0; i < sInfo.NumLights(); i++)
	{
		Vec3f lightDir;
		
		Color I = sInfo.GetLight(i)->Illuminate(sInfo, lightDir);
		lightDir.Normalize();

		Vec3f viewDir = sInfo.V();
		viewDir.Normalize();

		halfDir = lightDir + V;
		halfDir.Normalize();

		//=== DIFFUSE w/ ENERGY CONSERVATION ===//
		float diff = std::max(N.Dot(lightDir), 0.0f);
		Color diffuse = (diffColor * diff) / PI;

		//=== SPECULAR NORMALIZATION ZOOOOOOOOO ===//
		float spec = pow(std::max(N.Dot(halfDir), 0.0f), gloss) * ((gloss + 2.0) / (8.0 * PI * (1.0 - powf(2.0f, (-gloss / 2.0f) - 1.0f))));
		//float spec = pow(std::max(N.Dot(halfDir), 0.0f), gloss) * ((gloss + 2.0) / (4.0 * PI * (2.0 - exp(-gloss / 2.0))));
		//float spec = pow(std::max(N.Dot(halfDir), 0.0f), gloss);
		Color specular = specColor * spec;

		blinnPhong += I * (diffuse + specular);
	}

	float dist = BIGFLOAT;

	//============================== GLOBAL ILLUMINATION STUFF ==============================//
	const ShadeInfoExtend& extInfo = static_cast<const ShadeInfoExtend&>(sInfo);
	int giMax = GI_SAMPLES;

	//apple.EstimateIrradiance<50, PHOTONMAP_FILTER_QUADRATIC>(irrad, direction, 0.5f, sInfo.P(), sInfo.N(), 0.33f);
	//indirect = diffColor * irrad / PI;
	//orange.EstimateIrradiance<20, PHOTONMAP_FILTER_QUADRATIC>(irradCaustics, dirCaustics, 0.01f, sInfo.P(), sInfo.N(), 0.85f);
	//caustics = diffColor * irradCaustics / PI;
	
	
	//============== FUCKING IRRADIANCE STUFf ASDKLFHASDKLFHASDKLGHJASDKLGHASD ==============//
	// Makes sure that every random ray cast from global illumination only samples ONCE.
	if (extInfo.IsGI())
	{
		giMax = 1;
		apple.EstimateIrradiance<50, PHOTONMAP_FILTER_QUADRATIC>(irrad, direction, 0.05f, sInfo.P(), sInfo.N(), 0.33f);
		indirect = diffColor * irrad / PI;
	}

	else
	{
		orange.EstimateIrradiance<20, PHOTONMAP_FILTER_QUADRATIC>(irradCaustics, dirCaustics, 0.15f, sInfo.P(), sInfo.N(), 0.85f);
		caustics = diffColor * irradCaustics / PI;
	}


	// Put this into a loop in case I want to do a larger hemisphere sampling. Keep GI_SAMPLES at 1 for path tracing tho.
	if (extInfo.CurrentDepthGI() < GI_MAX)
	{
		extInfo.SetGI(true);
		extInfo.IncrementDepthGI();

		for (int n = 0; n < giMax; n++)
		{
			float rngX = sInfo.RandomFloat();
			float rngPhi = sInfo.RandomFloat() * 2.0f * PI;
			float hemiCos = Sqrt(1.0f - rngX);
			float hemiSin = Sqrt(rngX);
			//float hemiSin = sqrtf(1.0f - hemiCos * hemiCos);
			Vec3f orthX;
			Vec3f orthY;
			N.GetOrthonormals(orthX, orthY);
			Vec3f giSample = N * hemiCos + orthX * hemiSin * cosf(rngPhi) + orthY * hemiSin * sinf(rngPhi);
			giSample.Normalize();

			GI += sInfo.TraceSecondaryRay(giSample, dist);
		}

		extInfo.DecrementDepthGI();
		extInfo.SetGI(false);
	}

	GI /= giMax; //Normalize it so it doesn't blow up lol
	GI *= sInfo.Eval(diffuse);
	dist = BIGFLOAT;
	
	//=== FRESNEL ===//
	float cosThetaF = Clamp(fabs(N.Dot(V)), 0.0f, 1.0f);
	float F0 = pow(((1.0f - ior) / (1.0f + ior)), 2.0);
	float fresnel = F0 + (1.0f - F0) * pow(1.0f - cosThetaF, 5.0f);

	//=== REFRACTIONS PART ===//
	if (!refraction.GetValue().IsBlack())
	{
		float eta = 1.0f / ior;

		if (!sInfo.IsFront())
		{
			eta = 1 / eta;
		}

		Vec3f wo = sInfo.V();
		wo.Normalize();
		float oCos = wo.Dot(microfacetN);
		float tCosSq = 1.0f - (eta * eta) * (1.0f - (oCos * oCos));
		float tCos = Sqrt(Max(tCosSq, 0.0f));

		//=== ABSORPTION...FIX THIS ===//
		float absorbR = exp(-absorption.r * abs(sInfo.Depth() - sInfo.Depth()));
		float absorbG = exp(-absorption.g * abs(sInfo.Depth() - sInfo.Depth()));
		float absorbB = exp(-absorption.b * abs(sInfo.Depth() - sInfo.Depth()));
		Color absorb(absorbR, absorbG, absorbB);

		//=== ACTUAL REFRACTIONS ===//
		if (tCosSq >= 0.0f && sInfo.GetNode()->GetNodeObj() != nullptr)
		{
			Ray transRay;
			transRay.dir = (-eta * wo) - microfacetN * (tCos - eta * oCos);
			transRay.dir.Normalize();
			transRay.p = sInfo.P() + transRay.dir * BIAS;
			refractColor = sInfo.TraceSecondaryRay(transRay, dist) * refraction.GetValue() * absorb * (1.0 - fresnel);
		}

		//=== TOTAL INTERNAL REFLECTION - I think this is wrong???? ===//
		else
		{
			Ray internalRay;
			internalRay.dir = 2.0f * (sInfo.V().Dot(microfacetN)) * microfacetN - sInfo.V();
			internalRay.dir.Normalize();
			internalRay.p = sInfo.P() + internalRay.dir * BIAS;

			refractColor = sInfo.TraceSecondaryRay(internalRay, dist) * refraction.GetValue() * (1.0 - fresnel);
		}

		Ray internalRay;
		internalRay.dir = ReflectRay(sInfo, microfacetN);
		internalRay.p = sInfo.P() + internalRay.dir * BIAS;
		reflectColor = sInfo.TraceSecondaryRay(internalRay, dist) * absorb * refraction.GetValue() * fresnel;
		blinnPhong += refractColor + reflectColor;
	}

	//=== REFLECTIONS ===//
	if (!reflection.GetValue().IsBlack())
	{
		Ray reflectRay;
		reflectRay.dir = ReflectRay(sInfo, microfacetN);
		reflectRay.p = sInfo.P() + microfacetN * BIAS;
		reflectColor = sInfo.TraceSecondaryRay(reflectRay, dist) * reflection.GetValue();
		blinnPhong += reflectColor;
	}

	return indirect + caustics + blinnPhong + GI; //+ sInfo.Eval(emission) 
}

bool MtlBlinn::GenerateSample(SamplerInfo const& sInfo, Vec3f& dir, Info& si) const
{
	float gloss = sInfo.Eval(glossiness);
	float xRand = sInfo.RandomFloat();
	float phi = sInfo.RandomFloat() * 2.0f * PI;
	float cosTheta = powf(1.0f - xRand, 1.0f / (gloss + 1.0f));
	float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
	Vec3f u;
	Vec3f v;
	sInfo.N().GetOrthonormals(u, v);
	Vec3f microfacetN = sInfo.N() * cosTheta + u * sinTheta * cosf(phi) + v * sinTheta * sinf(phi);
	microfacetN.Normalize();

	Color d = sInfo.Eval(diffuse);
	Color s = sInfo.Eval(specular);
	Color r = sInfo.Eval(refraction);

	float rSum = d.r + s.r + r.r;
	float gSum = d.g + s.g + r.g;
	float bSum = d.b + s.b + r.b;

	float diffSum = d.r + d.g + d.b;
	float specSum = s.r + s.g + s.b;
	float refractSum = r.r + r.g + r.b;

	float totalSum = diffSum + specSum + refractSum;

	float Pr = Max(rSum, gSum, bSum);

	float Pd = (diffSum / totalSum) * Pr;
	float Ps = (specSum / totalSum) * Pr;
	float Pt = (refractSum / totalSum) * Pr;

	float totalP = Pd + Ps + Pt;
	Pd /= totalP;
	Ps /= totalP;
	Pt /= totalP;

	float rr = sInfo.RandomFloat();

	if (diffSum + specSum + refractSum <= 0.0f)
	{
		si.SetVoid();
		return false;
	}

	else if (rr <= Pd)
	{
		//hemi sample
		dir = CosWeightedHemisphereSampling(sInfo, sInfo.N());
		dir.Normalize();

		float cosi = std::max(EPSILON, sInfo.N().Dot(dir));

		si.lobe = DirSampler::DIFFUSE;
		si.mult = d / PI;
		si.prob = Pd * cosi / PI;
		return true;
	}

	else if (rr <= Pd + Ps)
	{
		//glossy fresnel
		float cosThetaH = std::max(0.0f, microfacetN.Dot(sInfo.V())); //F(wo, h)
		float F0 = pow(((1.0f - ior) / (1.0f + ior)), 2.0);
		float fresnelH = F0 + (1.0f - F0) * pow(1.0f - cosThetaH, 5.0f);

		float specNorm = pow(std::max(EPSILON, sInfo.N().Dot(microfacetN)), gloss) * ((gloss + 2.0) / (8.0 * PI * (1.0 - powf(2.0f, (-gloss / 2.0f) - 1.0f))));
		float hPDF = pow(std::max(EPSILON, sInfo.N().Dot(microfacetN)), gloss) * (gloss + 1.0f) / (2.0f * PI); //half vector PDF
		//float GGXndf = 

		//perfect reflection
		dir = ReflectRay(sInfo, microfacetN);
		dir.Normalize();
		
		float cosi = sInfo.N().Dot(dir);
		float coso = 1.0f / (sInfo.N().Dot(sInfo.V()));
		float J = 1.0f / (4.0f * sInfo.V().Dot(microfacetN));

		if (cosi < EPSILON || coso < EPSILON || J < EPSILON)
		{
			si.SetVoid();
			return false; 
		}

		si.lobe = DirSampler::SPECULAR;
		si.mult = s * specNorm / std::max(sInfo.V().Dot(microfacetN), EPSILON);
		si.prob = Ps * hPDF * J; 

		return true;
	}

	else if (rr <= Pd + Ps + Pt)
	{

		si.lobe = DirSampler::TRANSMISSION;
		si.mult = r;
		si.prob = Pt;

		if (sInfo.RandomFloat() < 0.02f)
		{
			si.mult = r * 2.0f;
			si.prob = 0.02f * 0.5f;
		}

		//refract
		Vec3f N = microfacetN;
		float etaT = this->IOR();

		const ShadeInfoExtend& extInfo = static_cast<const ShadeInfoExtend&>(sInfo);

		//WAVEY TIME
		float lambda = extInfo.GetPathLambda();
		float dispersionStrength = 0.5f; // 0.01 min - 0.05 max for rainbows
		float dispersiveIOR = GetDispersiveIOR(etaT, lambda, dispersionStrength);
		//=========

		float eta = 1.0f / dispersiveIOR;

		if (!sInfo.IsFront())
		{
			eta = 1 / eta;
			N = -N;
		}

		//fresnel
		Vec3f wo = sInfo.V();
		wo.Normalize();

		float cosThetaF = Clamp(fabs(sInfo.N().Dot(wo)), 0.0f, 1.0f);
		float F0 = pow(((1.0f - ior) / (1.0f + ior)), 2.0f);
		float fresnel = F0 + (1.0f - F0) * pow(1.0f - cosThetaF, 5.0f);

		float refractOrReflect = sInfo.RandomFloat();

		float oCos = N.Dot(wo);
		float tCosSq = 1.0f - (eta * eta) * (1.0f - (oCos * oCos));
		float tCos = Sqrt(tCosSq);
		

		//float Jt = fabs(tCos / powf((tCos + eta * oCos), 2.0f));
		float Jt = (eta * eta * fabs(oCos)) / fabs(tCos);
		float Js = 1.0f / (4.0f * sInfo.V().Dot(N));
		float hPDF = pow(std::max(EPSILON, sInfo.N().Dot(N)), gloss) * (gloss + 2.0f) / (2.0f * PI); //half vector PDF
		float specNorm = pow(std::max(sInfo.N().Dot(N), EPSILON), gloss) * ((gloss + 2.0) / (8.0 * PI * (1.0 - powf(2.0f, (-gloss / 2.0f) - 1.0f))));
		Color absorb = Color(1.0, 1.0, 1.0);

		if (refractOrReflect <= 1.0 - fresnel)
		{
			if (tCosSq >= 0.0f)
			{

				dir = eta * (-wo) + N * (eta * oCos - tCos);

				HitInfo absorbInfo;
				absorbInfo.z = std::numeric_limits<float>::max();
				absorbInfo.node = firstNode;

				Ray absorbRay;
				absorbRay.dir = dir;
				absorbRay.p = sInfo.P();

				NodeTraversal(firstNode, absorbRay , absorbInfo, 0);

				float distanceInMedium = fabs((sInfo.P() - absorbInfo.z).Length());
				absorb.r = exp(-absorption.r * distanceInMedium);
				absorb.g = exp(-absorption.g * distanceInMedium);
				absorb.b = exp(-absorption.b * distanceInMedium);


				si.mult = r * (1.0f - fresnel) * eta * eta * absorb;
				si.prob = Pt * (1.0f - fresnel);// * Jt;// * hPDF;
			}

			else
			{
				dir = ReflectRay(sInfo, N);
				si.mult = r;// *specNorm / std::max(dir.Dot(microfacetN), EPSILON);
				si.prob = Pt;// *hPDF* Js;

			}
		}

		else
		{
			dir = ReflectRay(sInfo, N);
			si.mult = r * fresnel;// *specNorm / std::max(dir.Dot(microfacetN), EPSILON);
			si.prob = Pt * fresnel;// *hPDF* Js;

		}

		dir.Normalize();
		return true;
	}

	else
	{
		si.SetVoid();
		return false;
	}
}

void MtlBlinn::GetSampleInfo(SamplerInfo const& sInfo, Vec3f const& dir, Info& si) const
{
	float gloss = sInfo.Eval(glossiness);

	Vec3f H = (sInfo.V() + dir);
	H.Normalize();

	Color d = sInfo.Eval(diffuse);
	Color s = sInfo.Eval(specular);
	Color r = sInfo.Eval(refraction);

	float rSum = d.r + s.r + r.r;
	float gSum = d.g + s.g + r.g;
	float bSum = d.b + s.b + r.b;

	float diffSum = d.r + d.g + d.b;
	float specSum = s.r + s.g + s.b;
	float refractSum = r.r + r.g + r.b;

	float totalSum = diffSum + specSum + refractSum;

	float Pr = Max(rSum, gSum, bSum);

	float Pd = (diffSum / totalSum) * Pr;
	float Ps = (specSum / totalSum) * Pr;
	float Pt = (refractSum / totalSum) * Pr;

	float totalP = Pd + Ps + Pt;
	Pd /= totalP;
	Ps /= totalP;
	Pt /= totalP;

	float rr = sInfo.RandomFloat();
	
	float cosLi = fabs(dir.Dot(sInfo.N()));

	float J = 1.0f / 4.0f * sInfo.V().Dot(H);
	float specNorm = pow(std::max(sInfo.N().Dot(H), 0.0f), gloss) * ((gloss + 2.0) / (8.0 * PI * (1.0 - powf(2.0f, (-gloss / 2.0f) - 1.0f))));
	float hPDF = pow(std::max(EPSILON, sInfo.N().Dot(H)), gloss) * (gloss + 2.0f) / (2.0f * PI); //half vector PDF

	float dProb = Pd * cosLi / PI; // cos(wi) / PI
	float sProb = Ps * hPDF * J;
	float tProb = Pt;

	si.prob = dProb + sProb + tProb;

	return;
}


Color MtlMicrofacet::Shade(ShadeInfo const& sInfo) const
{
	Color pbr;
	return pbr;
}

bool MtlMicrofacet::GenerateSample(SamplerInfo const& sInfo, Vec3f& dir, Info& si) const
{
	return false;
}

void MtlMicrofacet::GetSampleInfo(SamplerInfo const& sInfo, Vec3f const& dir, Info& si) const
{

	return;
}
