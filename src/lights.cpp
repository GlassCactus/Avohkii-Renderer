#pragma once
#include "lights.h"

#define PI 3.1415926f
#define BIAS 0.0002f
#define EPSILON 0.00001f
#define SPP_MIN 4
#define SPP_MAX 16

extern Node* firstNode;
void NodeTraversal(Node* root, Ray camToPixel, HitInfo& hInfo, int hitSide);

Color PointLight::Illuminate(ShadeInfo const& sInfo, Vec3f& dir) const
{
	float shadows = 0.0;

	//sampling a partial hemisphere
	Vec3f P = sInfo.P();
	Vec3f L = position - sInfo.P();
	dir = L.GetNormalized();
	float length = L.Length();
	float r2 = (Sqrt(length * length - size * size) * size) / length;
	Vec3f x;
	Vec3f y;
	dir.GetOrthonormals(x, y); 
	float t = sInfo.RandomFloat();
	float r = sInfo.RandomFloat();

	float shadowBreak = SPP_MAX;

	for (int n = 0; n < SPP_MAX; n++)
	{
		float rngT = (Halton(n, 2) + t) * 2.0f * PI;
		float rngR = Halton(n, 3) + r;

		if (rngR > 1)
		{
			rngR -= 1;
		}

		Vec2f inCircle = Vec2f(cos(rngT), sin(rngT));
		inCircle *= Sqrt(rngR) * r2;
		Vec3f d = (position + x * inCircle.x + y * inCircle.y) - sInfo.P();

		shadows += sInfo.TraceShadowRay(d, 1);

		if (n == SPP_MIN - 1)
		{
			if ((shadows/SPP_MIN) == 0.0f || (shadows/SPP_MIN) == 1.0f)
			{
				shadowBreak = SPP_MIN;
				break;
			}
		}
	}

	float att = 1.0f / (length * length);
	return att * intensity * (shadows/shadowBreak);
}

bool PointLight::IntersectRay(Ray const& ray, HitInfo& hInfo, int hitSide) const
{
	cyMatrix3f T = cyMatrix3f(Vec3f(size, 0.0f, 0.0f), Vec3f(0.0f, size, 0.0f), Vec3f(0.0f, 0.0f, size));
	cyMatrix3f invT = T.GetInverse();
	Ray transRay;
	transRay.p = invT * (ray.p - position);
	transRay.dir = invT * ray.dir;
	transRay.dir.Normalize();

	float A = transRay.dir.Dot(transRay.dir);
	float B = 2.0 * transRay.dir.Dot(transRay.p);
	float C = (transRay.p).Dot(transRay.p) - 1.0;
	float delta = B * B - 4.0 * A * C;

	if (delta > 0)
	{
		double t1 = (-B - Sqrt(delta)) / (2.0 * A);
		double t2 = (-B + Sqrt(delta)) / (2.0 * A);

		if (t1 >= BIAS && t2 >= BIAS && t1 < hInfo.z && t2 < hInfo.z)
		{
			hInfo.z = Min(t1, t2);
			hInfo.front = true;
			if (hitSide == 1)
			{
				return true;
			}
		}

		else if (t1 >= BIAS && t1 < hInfo.z)
		{
			hInfo.z = t1;
			hInfo.front = true;
			if (hitSide == 1)
			{
				return true;
			}
		}

		else if (t2 >= BIAS && t2 < hInfo.z)
		{
			hInfo.z = t2;
			hInfo.front = false;
		}

		else
		{
			return false;
		}

		hInfo.p = transRay.p + transRay.dir * hInfo.z;
		hInfo.N = hInfo.p;
		hInfo.N.Normalize();

		float u = (atan2(hInfo.p.y, hInfo.p.x) / (2.0f * PI)) + 0.5f;
		float v = (asin(hInfo.p.z) / (PI)) + 0.5f;
		hInfo.uvw = Vec3f(u, v, 0.0f);

		return true;
	}

	return false;

}

bool PointLight::GenerateSample(SamplerInfo const& sInfo, Vec3f& dir, Info& si) const 
{
	float xSphere = sInfo.RandomFloat();
	float phiSphere = sInfo.RandomFloat() * 2.0f * PI;
	float cosTheta = 1.0f - 2.0f * xSphere;
	float sinTheta = Sqrt(1.0f - (cosTheta * cosTheta));
	float length = (sInfo.P() - position).Length();
	float r2 = (Sqrt(length * length - size * size) * size) / length;

	dir = Vec3f(sinTheta * cosf(phiSphere), sinTheta * sinf(phiSphere), cosTheta);
	dir.Normalize();

	if (dir.Dot(sInfo.N()) > 0)
	{
		dir = -dir;
	}

	Vec3f lightN = dir;

	Vec3f lightOnSurface = position + size * dir;
	float area = 2.0f * PI * size * size;

	si.dist = (lightOnSurface - sInfo.P()).Length();
	dir = lightOnSurface - sInfo.P();
	dir.Normalize();
	si.mult = Radiance(sInfo);
	si.prob = (si.dist * si.dist) / (area * std::max(EPSILON, lightN.Dot(-dir)));
	
	//Occlusion Test
	Ray shadowRay;
	shadowRay.dir = dir;
	shadowRay.p = sInfo.P();
	HitInfo shadowHit;
	shadowHit.z = si.dist - EPSILON;
	shadowHit.node = firstNode;
	NodeTraversal(firstNode, shadowRay, shadowHit, 0);

	if (shadowHit.z < si.dist - EPSILON)
	{
		return false;
	}

	return true;
}

void PointLight::GetSampleInfo(SamplerInfo const& sInfo, Vec3f const& dir, Info& si) const
{
	Vec3f lightN = sInfo.P() - position;
	lightN.Normalize();
	Vec3f lightOnSurface = position + size * dir;
	si.dist = (lightOnSurface - sInfo.P()).Length();
	float area = 4.0f * PI * size * size;
	si.mult = Radiance(sInfo);
	si.prob = (si.dist * si.dist) / ((area) * fabs(lightN.Dot(-dir)));// This just ends up being 1;
	return;
}

void PointLight::RandomPhoton(RNG& rng, Ray& r, Color& c) const
{
	return;
}
