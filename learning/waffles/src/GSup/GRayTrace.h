/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GRAYTRACE_H__
#define __GRAYTRACE_H__

#include <stddef.h>
#include <math.h>
#include "G3D.h"
#include "../GClasses/GRand.h"
#include <vector>

namespace GClasses {

class GRayTraceLight;
class GRayTraceMaterial;
class GRayTraceObject;
class GRayTraceRay;
class GRayTraceScene;
class GRayTraceTriMesh;
class GRayTraceBoundingBoxBase;
class GNodeHashTable;
class GImage;

/// This class represents a color. It's more precise than GColor, but takes up more
/// memory. Note that the ray tracer ignores the alpha channel because the material
/// specifies a unique transmission color.
class GRayTraceColor
{
public:
	G3DReal a, r, g, b;

	GRayTraceColor()
	: a(1), r(0), g(0), b(0)
	{
	}

	GRayTraceColor(GRayTraceColor* pThat)
	{
		a = pThat->a;
		r = pThat->r;
		g = pThat->g;
		b = pThat->b;
	}

	GRayTraceColor(unsigned int c);

	GRayTraceColor(G3DReal alpha, G3DReal red, G3DReal green, G3DReal blue)
	{
		set(alpha, red, green, blue);
	}

	/// Save to a text-based format
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Load from a text-based format
	void fromTwt(GTwtNode* pNode);

	bool isBlack()
	{
		if(r == 0 && g == 0 && b == 0)
			return true;
		else
			return false;
	}

	void set(G3DReal alpha, G3DReal red, G3DReal green, G3DReal blue)
	{
		a = alpha;
		r = red;
		g = green;
		b = blue;
	}

	void set(unsigned int c);

	void copy(GRayTraceColor* pThat)
	{
		a = pThat->a;
		r = pThat->r;
		g = pThat->g;
		b = pThat->b;
	}

	void add(GRayTraceColor* pThat)
	{
		a = MAX(a, pThat->a);
		r += pThat->r;
		g += pThat->g;
		b += pThat->b;
	}

	void clip()
	{
		r = MIN((G3DReal)1, r);
		g = MIN((G3DReal)1, g);
		b = MIN((G3DReal)1, b);
	}

	void multiply(G3DReal mag)
	{
		r *= mag;
		g *= mag;
		b *= mag;
	}

	void multiply(GRayTraceColor* pThat)
	{
		a *= pThat->a;
		r *= pThat->r;
		g *= pThat->g;
		b *= pThat->b;
	}

	unsigned int color();

	void makeSliderColor(float f, GRayTraceColor* pDiffuseColor);
};




/// Represents the camera for a ray tracing scene
class GRayTraceCamera : public GCamera
{
protected:
	G3DReal m_focalDistance;
	G3DReal m_lensDiameter;
	int m_maxDepth; // max number of ray bounces

public:
	GRayTraceCamera(int width, int height) : GCamera(width, height)
	{
		m_focalDistance = 0;
		m_lensDiameter = 0;
		m_maxDepth = 4;
	}

	GRayTraceCamera(GTwtNode* pNode);

	virtual ~GRayTraceCamera()
	{
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc);
	void setMaxDepth(int val) { m_maxDepth = val; }
	void setFocalDistance(G3DReal val) { m_focalDistance = val; }
	void setLensDiameter(G3DReal val) { m_lensDiameter = val; }
	inline int maxDepth() { return m_maxDepth; }
	inline G3DReal focalDistance() { return m_focalDistance; }
	inline G3DReal lensDiameter() { return m_lensDiameter; }
};



/// Represents a scene that you can ray-trace
class GRayTraceScene
{
public:
	enum RenderMode
	{
		FAST_RAY_TRACE,
		QUALITY_RAY_TRACE,
		PATH_TRACE,
	};

protected:
	GRayTraceColor m_backgroundColor;
	GRayTraceColor m_ambientLight;
	std::vector<GRayTraceMaterial*> m_materials;
	std::vector<GRayTraceTriMesh*> m_meshes; // meshes reference a material
	std::vector<GRayTraceObject*> m_objects; // triangles reference a mesh, spheres reference a material
	std::vector<GRayTraceLight*> m_lights; // area lights reference an object
	GRayTraceCamera* m_pCamera;
	GRayTraceBoundingBoxBase* m_pBoundingBoxTree;
	G3DReal m_toneMappingConstant;

	/// Rendering values
	RenderMode m_eMode;
	GImage* m_pImage;
	G3DReal* m_pDistanceMap;
	int m_nY;
	G3DVector m_pixSide;
	G3DVector m_pixDX;
	G3DVector m_pixDY;
	GRand* m_pRand;

public:
	GRayTraceScene(GRand* pRand);
	GRayTraceScene(GTwtNode* pNode, GRand* pRand);
	~GRayTraceScene();

	/// Save to a text-based format
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Deletes all the materials and objects. (Leaves the lights
	/// and camera as they are.)
	void flushObjects();

	/// Swap all matrials, meshes and objects with the other scene. (Leaves
	/// the camera and lights as they are.)
	void swapObjects(GRayTraceScene* pOther);

	/// Specify whether to emphasize quality or speed
	void setRenderMode(RenderMode eMode) { m_eMode = eMode; }

	/// This method calls RenderBegin, then calls RenderLine until
	/// the whole image has been rendered.
	void render();

	/// Call this before calling RenderLine(). It resets the image and
	/// computes values necessary for rendering.
	void renderBegin();

	/// Call this to render a singe horizontal line of the image. Returns true if there's
	/// still more rendering to do. Returns false if it's done. You must call RenderBegin()
	/// once before you start calling this method.
	bool renderLine();

	/// This draws a wire frame of the scene. This is a fast way to ensure your
	/// camera is looking where you think it is looking before you perform a long render.
	void drawWireFrame();

	/// This calls RenderBegine and then renders a single pixel. It's not efficient
	/// to call this method for every pixel. The only purpose for this method is to
	/// make debugging the ray tracer easier. Pick a pixel that isn't rendered the way
	/// you want and step through the ray tracing process to see why.
	unsigned int renderSinglePixel(int x, int y);

	/// Returns the rendered image (or partially rendered image). Returns NULL if Render()
	/// or RenderBegin() has not been called yet.
	GImage* image() { return m_pImage; }

	GImage* releaseImage()
	{
		GImage* pImage = m_pImage;
		m_pImage = NULL;
		return pImage;
	}

	void setBackgroundColor(G3DReal a, G3DReal r, G3DReal g, G3DReal b)
	{
		m_backgroundColor.set(a, r, g, b);
	}

	void setAmbientLight(G3DReal r, G3DReal g, G3DReal b)
	{
		m_ambientLight.set(1, r, g, b);
	}

	GRayTraceColor* ambientLight() { return &m_ambientLight; }

	GRayTraceCamera* camera()
	{
		return m_pCamera;
	}

	GRand* rand() { return m_pRand; }
	void activateDistanceMap();
	G3DReal* distanceMap() { return m_pDistanceMap; }
	void addMaterial(GRayTraceMaterial* pMaterial);
	void addMesh(GRayTraceTriMesh* pMesh);
	void addObject(GRayTraceObject* pObject);
	void addLight(GRayTraceLight* pLight);
	size_t materialCount();
	size_t meshCount();
	size_t objectCount();
	size_t lightCount();
	GRayTraceMaterial* material(size_t n);
	GRayTraceTriMesh* mesh(size_t n);
	GRayTraceObject* object(size_t n);
	GRayTraceLight* light(size_t n);
	std::vector<GRayTraceMaterial*>& materials() { return m_materials; }
	GRayTraceColor* backgroundColor() { return &m_backgroundColor; }
	GRayTraceBoundingBoxBase* boundingBoxTree() { return m_pBoundingBoxTree; }
	void setToneMappingConstant(G3DReal c) { m_toneMappingConstant = c; }
	size_t materialIndex(GRayTraceMaterial* pMaterial);
	size_t meshIndex(GRayTraceTriMesh* pMesh);
	size_t objectIndex(GRayTraceObject* pObj);

	unsigned int renderPixel(GRayTraceRay* pRay, G3DVector* pScreenPoint, G3DReal* pDistance);
	unsigned int renderPixelAntiAliassed(GRayTraceRay* pRay, G3DVector* pScreenPoint, G3DReal* pDistance);
	unsigned int renderPixelPathTrace(GRayTraceRay* pRay, G3DVector* pScreenPoint);
};



/// Represents a source of light in a ray-tracing scene
class GRayTraceLight
{
public:
	enum LightType
	{
		Directional = 0,
		Point,
		Area,
	};

protected:
	GRayTraceColor m_color;

public:
	GRayTraceLight(G3DReal r, G3DReal g, G3DReal b);
	GRayTraceLight(GTwtNode* pNode);
	virtual ~GRayTraceLight();

	static GRayTraceLight* fromTwt(GTwtNode* pNode, GRayTraceScene* pScene);
	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene) = 0;
	virtual LightType lightType() = 0;
	virtual void colorContribution(GRayTraceScene* pScene, GRayTraceRay* pRay, GRayTraceMaterial* pMaterial, bool bSpecular) = 0;

protected:
	GTwtNode* baseTwtNode(GTwtDoc* pDoc);
};




/// Represents directional light in a ray-tracing scene
class GRayTraceDirectionalLight : public GRayTraceLight
{
protected:
	G3DVector m_direction;
	G3DReal m_jitter;

public:
	/// Specify the direction to the light, not the direction the light travels
	GRayTraceDirectionalLight(G3DReal dx, G3DReal dy, G3DReal dz, G3DReal r, G3DReal g, G3DReal b, G3DReal jitter);
	GRayTraceDirectionalLight(GTwtNode* pNode);
	virtual ~GRayTraceDirectionalLight();

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene);
	virtual LightType lightType() { return Directional; }
	virtual void colorContribution(GRayTraceScene* pScene, GRayTraceRay* pRay, GRayTraceMaterial* pMaterial, bool bSpecular);
};



/// Represents a point light in a ray-tracing scene
class GRayTracePointLight : public GRayTraceLight
{
protected:
	G3DVector m_position;
	G3DReal m_jitter;

public:
	GRayTracePointLight(G3DReal x, G3DReal y, G3DReal z, G3DReal r, G3DReal g, G3DReal b, G3DReal jitter);
	GRayTracePointLight(GTwtNode* pNode);
	virtual ~GRayTracePointLight();

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene);
	virtual LightType lightType() { return Point; }
	virtual void colorContribution(GRayTraceScene* pScene, GRayTraceRay* pRay, GRayTraceMaterial* pMaterial, bool bSpecular);
};



/// Represents a light source with area
class GRayTraceAreaLight : public GRayTraceLight
{
protected:
	GRayTraceObject* m_pObject;

public:
	GRayTraceAreaLight(GRayTraceObject* pObject, G3DReal r, G3DReal g, G3DReal b);
	GRayTraceAreaLight(GTwtNode* pNode, GRayTraceScene* pScene);
	virtual ~GRayTraceAreaLight();

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene);
	virtual LightType lightType() { return Area; }
	virtual void colorContribution(GRayTraceScene* pScene, GRayTraceRay* pRay, GRayTraceMaterial* pMaterial, bool bSpecular);
};



class GRayTraceMaterial
{
public:
	enum MaterialType
	{
		Physical,
		Image,
		Etherial,
	};

	enum ColorType
	{
		Diffuse = 0, // if the ray can see the light from the point of intersection
		Specular, // if the ray reflects into a light
		Reflective, // colors reflected by the object
		Transmissive, // colors that pass through the object
		Ambient, // light that is always present in every direction
		Emissive, // for object lights or flourescent objects
		Color_Type_Count,
	};

	GRayTraceMaterial();
	virtual ~GRayTraceMaterial();

	static GRayTraceMaterial* fromTwt(GTwtNode* pNode);

	virtual GTwtNode* toTwt(GTwtDoc* pDoc) = 0;
	virtual MaterialType materialType() = 0;
	virtual GRayTraceColor* color(ColorType eType, GRayTraceRay* pRay) = 0;
	virtual G3DReal indexOfRefraction() = 0;
	virtual G3DReal specularExponent() = 0;
	virtual G3DReal glossiness() = 0; // how much to jitter reflected rays
	virtual G3DReal cloudiness() = 0; // how much to jitter transmitted rays
	virtual bool isSame(GRayTraceMaterial* pThat) = 0;
	virtual GRayTraceMaterial* copy() = 0;

	void computeColor(GRayTraceScene* pScene, GRayTraceRay* pRay, bool bAmbient, bool bSpecular);
};


/// Represents the material of which an object is made in a ray-tracing scene
class GRayTracePhysicalMaterial : public GRayTraceMaterial
{
protected:
	GRayTraceColor m_colors[Color_Type_Count];
	G3DReal m_indexOfRefraction;
	G3DReal m_specularExponent;
	G3DReal m_glossiness;
	G3DReal m_cloudiness;

public:
	GRayTracePhysicalMaterial();
	GRayTracePhysicalMaterial(GTwtNode* pNode);
	virtual ~GRayTracePhysicalMaterial();

	virtual GTwtNode* toTwt(GTwtDoc* pDoc);
	virtual MaterialType materialType() { return Physical; }

	/// Ignores pRay and returns the color of the specified type. (pRay is
	/// used by image-texture materials to determine which pixel applies.)
	virtual GRayTraceColor* color(ColorType eType, GRayTraceRay* pRay);

	virtual G3DReal indexOfRefraction() { return m_indexOfRefraction; }
	virtual G3DReal specularExponent() { return m_specularExponent; }
	virtual G3DReal glossiness() { return m_glossiness; }
	virtual G3DReal cloudiness() { return m_cloudiness; }
	virtual bool isSame(GRayTraceMaterial* pThat);
	virtual GRayTraceMaterial* copy();

	void setColor(ColorType eType, G3DReal r, G3DReal g, G3DReal b);
	void setColor(ColorType eType, GRayTraceColor* pCol);
	void setIndexOfRefraction(G3DReal val) { m_indexOfRefraction = val; }
	void setSpecularExponent(G3DReal val) { m_specularExponent = val; }
	void setGlossiness(G3DReal val) { m_glossiness = val; }
	void setCloudiness(G3DReal val) { m_cloudiness = val; }
};


class GRayTraceImageTexture : public GRayTraceMaterial
{
protected:
	GImage* m_pTextureImage;
	bool m_bDeleteTextureImage;
	GRayTraceColor m_col;

public:
	GRayTraceImageTexture();
	GRayTraceImageTexture(GTwtNode* pNode);
	virtual ~GRayTraceImageTexture();

	virtual GTwtNode* toTwt(GTwtDoc* pDoc);
	virtual MaterialType materialType() { return Image; }
	virtual GRayTraceColor* color(ColorType eType, GRayTraceRay* pRay);
	virtual G3DReal indexOfRefraction() { return 1; }
	virtual G3DReal specularExponent() { return 1; }
	virtual G3DReal glossiness() { return 0; }
	virtual G3DReal cloudiness() { return 0; }
	virtual bool isSame(GRayTraceMaterial* pThat);
	virtual GRayTraceMaterial* copy();

	void setTextureImage(GImage* pImage, bool bDeleteImage);
	GImage* textureImage() { return m_pTextureImage; }
};




#define MAX_OBJECTS_PER_BOUNDING_BOX 5

/// A class used for making ray-tracing faster
class GRayTraceBoundingBoxBase
{
public:
	G3DVector m_min;
	G3DVector m_max;

public:
	GRayTraceBoundingBoxBase() {}
	virtual ~GRayTraceBoundingBoxBase() {}

	static GRayTraceBoundingBoxBase* makeBoundingBoxTree(GRayTraceScene* pScene);
	virtual bool isLeaf() = 0;
	virtual GRayTraceObject* closestIntersection(G3DVector* pRayOrigin, G3DVector* pDirectionVector, G3DReal* pOutDistance) = 0;

protected:
	static GRayTraceBoundingBoxBase* BuildTree(std::vector<GRayTraceObject*>& objects);
	bool DoesRayHitBox(G3DVector* pRayOrigin, G3DVector* pDirectionVector);
};


/// A class used for making ray-tracing faster
class GRayTraceBoundingBoxInterior : public GRayTraceBoundingBoxBase
{
protected:
	GRayTraceBoundingBoxBase* m_pLesser;
	GRayTraceBoundingBoxBase* m_pGreater;

public:
	GRayTraceBoundingBoxInterior(GRayTraceBoundingBoxBase* pLesser, GRayTraceBoundingBoxBase* pGreater)
	: GRayTraceBoundingBoxBase()
	{
		m_pLesser = pLesser;
		m_pGreater = pGreater;
		m_min.m_vals[0] = MIN(pLesser->m_min.m_vals[0], pGreater->m_min.m_vals[0]);
		m_min.m_vals[1] = MIN(pLesser->m_min.m_vals[1], pGreater->m_min.m_vals[1]);
		m_min.m_vals[2] = MIN(pLesser->m_min.m_vals[2], pGreater->m_min.m_vals[2]);
		m_max.m_vals[0] = MAX(pLesser->m_max.m_vals[0], pGreater->m_max.m_vals[0]);
		m_max.m_vals[1] = MAX(pLesser->m_max.m_vals[1], pGreater->m_max.m_vals[1]);
		m_max.m_vals[2] = MAX(pLesser->m_max.m_vals[2], pGreater->m_max.m_vals[2]);
	}

	virtual ~GRayTraceBoundingBoxInterior()
	{
		delete(m_pLesser);
		delete(m_pGreater);
	}

	virtual bool isLeaf() { return false; };
	virtual GRayTraceObject* closestIntersection(G3DVector* pRayOrigin, G3DVector* pDirectionVector, G3DReal* pOutDistance);
};


/// A class used for making ray-tracing faster
class GRayTraceBoundingBoxLeaf : public GRayTraceBoundingBoxBase
{
protected:
	int m_nObjectCount;
	GRayTraceObject** m_pObjects;

public:
	GRayTraceBoundingBoxLeaf(std::vector<GRayTraceObject*>& objects);
	virtual ~GRayTraceBoundingBoxLeaf();

	virtual bool isLeaf() { return true; }
	virtual GRayTraceObject* closestIntersection(G3DVector* pRayOrigin, G3DVector* pDirectionVector, G3DReal* pOutDistance);
};



/// Represents a triangle mesh in a ray-tracing scene
class GRayTraceTriMesh
{
protected:
	GRayTraceMaterial* m_pMaterial;
	int m_nPoints;
	G3DVector* m_pPoints;
	int m_nTriangles;
	int* m_pTriangles;
	G3DVector* m_pNormals;
	G3DReal* m_pTextureCoords;
	bool m_bCulling;

public:
	GRayTraceTriMesh(GRayTraceMaterial* pMaterial, int nPoints, int nTriangles, int nNormals, int nTextureCoords);
	GRayTraceTriMesh(GTwtNode* pNode, GRayTraceScene* pScene);
	~GRayTraceTriMesh();

	GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene);

	static GRayTraceTriMesh* makeCylinder(GRayTraceMaterial* pMaterial, G3DVector* pCenter1, G3DVector* pCenter2, G3DReal radius, int nSides, bool bEndCaps);

	/// The vertices must go around the surface. Don't cut across corners.
	static GRayTraceTriMesh* makeQuadSurface(GRayTraceMaterial* pMaterial, G3DVector* p1, G3DVector* p2, G3DVector* p3, G3DVector* p4);

	static GRayTraceTriMesh* makeSingleTriangle(GRayTraceMaterial* pMaterial, G3DVector* p1, G3DVector* p2, G3DVector* p3);

	G3DReal rayDistanceToTriangle(int nTriangle, G3DVector* pRayOrigin, G3DVector* pRayDirection);
	void normalVector(GRayTraceRay* pRay, int nIndex);
	bool isCulled() { return m_bCulling || m_pNormals; }
	void activateCulling() { m_bCulling = true; }
	void setPoint(int nIndex, const G3DVector* pPoint);
	void setTriangle(int nIndex, int v1, int v2, int v3);
	void setNormal(int nIndex, G3DVector* pNormal);
	void setTextureCoord(int nIndex, G3DReal x, G3DReal y);
	void center(G3DVector* pOutPoint, int nIndex);
	int triangleCount() { return m_nTriangles; }
	void triangle(int index, int* v1, int* v2, int* v3);
	GRayTraceMaterial* material() { return m_pMaterial; }
	void adjustBoundingBox(int nIndex, G3DVector* pMin, G3DVector* pMax);
	G3DVector* vertex(int nIndex, int nVertex);

	/// Automatically compute phong normals at the vertices to make the object appear smooth
	void computePhongNormals();

protected:
	bool isPointWithinPlanarPolygon(G3DVector* pPoint, G3DVector** ppVertices, int nVertices);
	//void ComputeTriangleNeighborSides(GIntQueue* pQ, GNodeHashTable* pEdges, int* pFaces, int i);
};




/// An object in a ray-tracing scene
class GRayTraceObject
{
public:
	enum ObjectType
	{
		Sphere,
		Triangle,
	};

	GRayTraceObject()
	{
	}

	virtual ~GRayTraceObject()
	{
	}

	static GRayTraceObject* fromTwt(GTwtNode* pNode, GRayTraceScene* pScene);
	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene) = 0;
	virtual ObjectType type() = 0;
	virtual GRayTraceMaterial* material() = 0;
	virtual G3DReal rayDistance(G3DVector* pRayOrigin, G3DVector* pRayDirection) = 0;
	virtual void normalVector(GRayTraceRay* pRay) = 0;
	virtual bool isCulled() = 0;
	virtual void center(G3DVector* pOutPoint) = 0;
	virtual void adjustBoundingBox(G3DVector* pMin, G3DVector* pMax) = 0;
	virtual void drawWireFrame(GCamera* pCamera, GImage* pImage) = 0;
};



/// A sphere in a ray-tracing scene
class GRayTraceSphere : public GRayTraceObject
{
protected:
	GRayTraceMaterial* m_pMaterial;
	G3DVector m_center;
	G3DReal m_radius;

public:
	GRayTraceSphere(GRayTraceMaterial* pMaterial, G3DReal x, G3DReal y, G3DReal z, G3DReal radius)
	: GRayTraceObject(), m_pMaterial(pMaterial), m_center(x, y, z), m_radius(radius)
	{
	}

	GRayTraceSphere(GTwtNode* pNode, GRayTraceScene* pScene);

	virtual ~GRayTraceSphere()
	{
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene);
	virtual ObjectType type() { return Sphere; }
	virtual GRayTraceMaterial* material() { return m_pMaterial; }
	virtual G3DReal rayDistance(G3DVector* pRayOrigin, G3DVector* pRayDirection);
	virtual void normalVector(GRayTraceRay* pRay);
	virtual bool isCulled() { return false; }
	virtual void center(G3DVector* pOutPoint);
	virtual void adjustBoundingBox(G3DVector* pMin, G3DVector* pMax);
	virtual void drawWireFrame(GCamera* pCamera, GImage* pImage);
	G3DVector* center() { return &m_center; }
	G3DReal radius() { return m_radius; }
};





/// A single triangle in a ray-tracing scene
class GRayTraceTriangle : public GRayTraceObject
{
protected:
	GRayTraceTriMesh* m_pMesh;
	int m_nIndex;

public:
	GRayTraceTriangle(GRayTraceTriMesh* pMesh, int nIndex)
	: GRayTraceObject(), m_pMesh(pMesh), m_nIndex(nIndex)
	{
	}

	GRayTraceTriangle(GTwtNode* pNode, GRayTraceScene* pScene);

	virtual ~GRayTraceTriangle()
	{
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, GRayTraceScene* pScene);
	virtual ObjectType type() { return Triangle; }
	virtual G3DReal rayDistance(G3DVector* pRayOrigin, G3DVector* pRayDirection);
	virtual void normalVector(GRayTraceRay* pRay);
	virtual bool isCulled() { return m_pMesh->isCulled(); }
	virtual GRayTraceMaterial* material() { return m_pMesh->material(); }
	virtual void center(G3DVector* pOutPoint);
	virtual void adjustBoundingBox(G3DVector* pMin, G3DVector* pMax);
	virtual void drawWireFrame(GCamera* pCamera, GImage* pImage);
	G3DVector* vertex(int nVertex);
};


} // namespace GClasses

#endif // __GRAYTRACE_H__
