#include <iostream>
#include <filesystem>
#include <thread>

#include "renderer.h"
#include "scene.h"
#include "objects.h"
#include "tinyxml2.h"
#include "materials.h"
#include "lights.h"

#define PI 3.1415926f
#define BIAS 0.0002f

//const char* filename = "C:\\Users\\aznbo\\Desktop\\Av-Renderer\\src\\chainsaw.xml";
const char* filename = "C:\\Users\\aznbo\\Desktop\\Av-Renderer\\src\\ajax.xml";


bool Plane::IntersectRay(Ray const& ray, HitInfo& hInfo, int hitSide) const
{
	Vec3f planeN = Vec3f(0.0f, 0.0f, 1.0f);
	float dotNormal = ray.dir.Dot(planeN);

	if (fabs(dotNormal) <= BIAS)
	{			
		return false;
	}

	float t = -ray.p.Dot(planeN) / dotNormal;

	if (t > BIAS && t < hInfo.z)
	{
		Vec3f tempHit = ray.p + ray.dir * t;

		if (tempHit.x > 1 || tempHit.x < -1 || tempHit.y > 1 || tempHit.y < -1)
		{
			return false;
		}

		hInfo.z = t;
		hInfo.p = tempHit;

		float u = 0.5f * (hInfo.p.x + 1.0f);
		float v = 0.5f * (hInfo.p.y + 1.0f);
		hInfo.uvw = Vec3f(u, v, 0.0f);

		hInfo.front = (dotNormal < 0);
		hInfo.N = planeN; // hInfo.front ? planeN : -planeN; this is taken care of in materials.xpp
		hInfo.N.Normalize();

		return true;
	}

	return false;
}

//========== AABB INTERSECTION ==========//
bool Box::IntersectRay(Ray const& r, float t_max) const
{
	Vec3f P = r.p;
	Vec3f D = r.dir;

	float txMax = (pmax[0] - P[0]) / D[0];
	float txMin = (pmin[0] - P[0]) / D[0];
	if (txMin > txMax)
	{
		std::swap(txMax, txMin);
	}

	float tyMax = (pmax[1] - P[1]) / D[1];
	float tyMin = (pmin[1] - P[1]) / D[1];
	if (tyMin > tyMax)
	{
		std::swap(tyMax, tyMin);
	}

	float tzMax = (pmax[2] - P[2]) / D[2];
	float tzMin = (pmin[2] - P[2]) / D[2];
	if (tzMin > tzMax)
	{
		std::swap(tzMax, tzMin);
	}

	float tMax = Min(Min(txMax, tyMax), tzMax);
	float tMin = Max(Max(txMin, tyMin), tzMin);

	if (tMin > tMax || tMax < 0)
	{
		return false;
	}

	return true;
}

bool TriObj::TraceBVHNode(Ray const& ray, HitInfo& hInfo, int hitSide, unsigned int nodeID) const
{
	bool hit = false;

	if (!bvh.IsLeafNode(nodeID))
	{
		unsigned int left = bvh.GetFirstChildNode(nodeID);
		unsigned int right = bvh.GetSecondChildNode(nodeID);
		Box leftBox = bvh.GetNodeBounds(left);
		Box rightBox = bvh.GetNodeBounds(right);

		if (leftBox.IntersectRay(ray, 0))
		{
			hit |= TraceBVHNode(ray, hInfo, hitSide, left);
		}

		if (rightBox.IntersectRay(ray, 0))
		{
			hit |= TraceBVHNode(ray, hInfo, hitSide, right);
		}

		return hit;
	}

	int triangleCount = bvh.GetNodeElementCount(nodeID);
	unsigned int const* triangleInBox = bvh.GetNodeElements(nodeID);

	for (int i = 0; i < triangleCount; i++)
	{
		hit |= IntersectTriangle(ray, hInfo, hitSide, triangleInBox[i]);
	}

	return hit;
}

bool TriObj::IntersectTriangle(Ray const& ray, HitInfo& hInfo, int hitSide, unsigned int faceID) const
{
	float epsilon = 0.00002f;
	Vec3f D = ray.dir;
	Vec3f P = ray.p;

	TriFace f = F(faceID);
	Vec3f A = V(f.v[0]);
	Vec3f B = V(f.v[1]);
	Vec3f C = V(f.v[2]);

	Vec3f T = P - A;
	Vec3f E1 = B - A;
	Vec3f E2 = C - A;
	Vec3f Ecross = E1.Cross(E2);

	float detM = (-D).Dot(Ecross);

	if (fabs(detM) < epsilon)
	{
		return false;
	}

	float Xdet = T.Dot(Ecross);
	float Ydet = (-D).Dot(T.Cross(E2));
	float Zdet = (-D).Dot(E1.Cross(T));

	float t = Xdet / detM;
	float u = Ydet / detM;
	float v = Zdet / detM;

	if (t < -epsilon || t > hInfo.z)
	{
		return false;
	}

	if (u < -epsilon || v < -epsilon || u + v > 1.0f + epsilon)
	{
		return false;
	}

	if (t > epsilon && t <= hInfo.z)
	{
		hInfo.z = t;
		hInfo.p = P + D * hInfo.z;
		hInfo.N = GetNormal(faceID, Vec3f(1.0 - u - v, u, v));
		hInfo.N.Normalize();

		if (this->HasTextureVertices())
		{
			Vec3f texCoordA = this->GetTexCoord(faceID, Vec3f(1.0f, 0.0f, 0.0f));
			Vec3f texCoordB = this->GetTexCoord(faceID, Vec3f(0.0f, 1.0, 0.0f));
			Vec3f texCoordC = this->GetTexCoord(faceID, Vec3f(0.0f, 0.0f, 1.0f));

			hInfo.uvw = texCoordA * (1.0f - u - v) + texCoordB * u + texCoordC * v;
		}

		hInfo.front = (D.Dot(hInfo.N) <= 0);

		if (hitSide == 1) //I'm using hitSide for shadows. If 1, then return on literally ANY hit.
		{
			return true;
		}

		return true;
	}


	if (hInfo.z >= std::numeric_limits<float>::max() - epsilon)
	{
		return false;
	}

	return false;
}

//=============== TROMBONES ===============//
bool TriObj::IntersectRay(Ray const& ray, HitInfo& hInfo, int hitSide) const
{
	//hInfo.z = std::numeric_limits<float>::max();//why the fuck do I have this???
	return TraceBVHNode(ray, hInfo, hitSide, bvh.GetRootNodeID());
}


//AABB Just leaving this here in case I need a reminder on how it's calculated lol.
/*
if (AABB.IsEmpty())
{
	std::cout << "prepping ingredients" << std::endl;

	for (int i = 0; i < nf; i++) //goes through all faces
	{
		TriMesh::TriFace b = F(i);
		Vec3f aBox = V(b.v[0]);
		Vec3f bBox = V(b.v[1]);
		Vec3f cBox = V(b.v[2]);

		for (int j = 0; j < 3; j++) //checks through all vertices and find the max and min coordinates
		{
			float tempMax = Max(Max(aBox[j], bBox[j]), cBox[j]);
			if (tempMax >= AABB.pmax[j])
			{
				AABB.pmax[j] = tempMax;
			}

			float tempMin = Min(Min(aBox[j], bBox[j]), cBox[j]);
			if (tempMin <= AABB.pmin[j])
			{
				AABB.pmin[j] = tempMin;
			}
		}
	}

	std::cout << "done cooking" << std::endl;
}
*/

bool Sphere::IntersectRay(Ray const& ray, HitInfo& hInfo, int hitSide) const
{
	float A = ray.dir.Dot(ray.dir);
	float B = 2.0 * ray.dir.Dot(ray.p);
	float C = (ray.p).Dot(ray.p) - 1.0;
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

		hInfo.p = ray.p + ray.dir * hInfo.z;
		hInfo.N = hInfo.p;
		hInfo.N.Normalize();

		float u = (atan2(hInfo.p.y, hInfo.p.x) / (2.0f * PI)) + 0.5f;
		float v = (asin(hInfo.p.z) / (PI)) + 0.5f;
		hInfo.uvw = Vec3f(u, v, 0.0f);

		return true;
	}

	return false;

}


//This will always be checking the select node's CHILD. IF the node doesn't have a child it was return but that's ok b/c it will be checked in the prev loop.
//The RAYS will all be FROM the camera TO the pixel centers in case you ever forget
void NodeTraversal(Node* root, Ray camToPixel, HitInfo& hInfo, int hitSide)
{
	if (root->GetNumChild() <= 0)
	{
		return;
	}

	for (int j = 0; j < root->GetNumChild(); j++)
	{
		Node* child = root->GetChild(j);
		Ray rLocal = child->ToNodeCoords(camToPixel);
		NodeTraversal(child, rLocal, hInfo, hitSide);

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


int main()
{
	RayTracer renderer;
	renderer.LoadScene(filename);
	ShowViewport(&renderer, false);
	return 0;
}