#ifndef _SUBDIV_2D_H_
#define _SUBDIV_2D_H_

#include "amcomdef.h"
#include <vector>

typedef struct
{
	MFloat pt1x;
	MFloat pt1y;
	MFloat pt2x;
	MFloat pt2y;
	MFloat pt3x;
	MFloat pt3y;
} VEC6F;

typedef struct
{
    MFloat x;
    MFloat y;
}MPOINTF;

/** Subdiv2D point location cases */
enum {
	PTLOC_ERROR = -2, //!< Point location error
	PTLOC_OUTSIDE_RECT = -1, //!< Point outside the subdivision bounding rect
	PTLOC_INSIDE = 0, //!< Point inside some facet
	PTLOC_VERTEX = 1, //!< Point coincides with one of the subdivision vertices
	PTLOC_ON_EDGE = 2  //!< Point on some edge
};

/** Subdiv2D edge type navigation (see: getEdge()) */
enum {
	NEXT_AROUND_ORG = 0x00,
	NEXT_AROUND_DST = 0x22,
	PREV_AROUND_ORG = 0x11,
	PREV_AROUND_DST = 0x33,
	NEXT_AROUND_LEFT = 0x13,
	NEXT_AROUND_RIGHT = 0x31,
	PREV_AROUND_LEFT = 0x20,
	PREV_AROUND_RIGHT = 0x02
};


class Subdiv2D
{
public:

	Subdiv2D();
	~Subdiv2D();

	void initDelaunay(MRECTF rect);
	int insert(MPOINTF pt);
	void getTriangleList(std::vector<VEC6F>& triangleList) const;

	int locate(MPOINTF pt, int& _edge, int& _vertex);
	int rotateEdge(int edge, int rotate) const;
	int nextEdge(int edge) const;
	int symEdge(int edge) const;
	int getEdge(int edge, int nextEdgeType) const;
	MPOINTF getVertex(int vertex, int* firstEdge = 0) const;
	int edgeOrg(int edge, MPOINTF* orgpt = 0) const;
	int edgeDst(int edge, MPOINTF* dstpt = 0) const;

private:
	int newEdge();
	void deleteEdge(int edge);
	int newPoint(MPOINTF pt, bool isvirtual, int firstEdge = 0);
	void deletePoint(int vtx);
	void setEdgePoints(int edge, int orgPt, int dstPt);
	void splice(int edgeA, int edgeB);
	int connectEdges(int edgeA, int edgeB);
	void swapEdges(int edge);
	int isRightOf(MPOINTF pt, int edge) const;
	//void calcVoronoi();
	//void clearVoronoi();
	//void checkSubdiv() const;
	struct Vertex
	{
		Vertex();
		Vertex(MPOINTF pt, bool _isvirtual, int _firstEdge = 0);
		bool isvirtual() const;
		bool isfree() const;

		int firstEdge;
		int type;
		MPOINTF pt;
	};

	struct QuadEdge
	{
		QuadEdge();
		QuadEdge(int edgeidx);
		bool isfree() const;

		int next[4];
		int pt[4];
	};

private:
	//! All of the vertices
	std::vector<Vertex> vtx;
	//! All of the edges
	std::vector<QuadEdge> qedges;
	int freeQEdge;
	int freePoint;
	bool validGeometry;

	int recentEdge;
	//! Top left corner of the bounding rect
	MPOINTF topLeft;
	//! Bottom right corner of the bounding rect
	MPOINTF bottomRight;

};

#endif // !_SUBDIV_2D_H_

