#include "delTriangulation.h"
#include <cstddef>
DelTriangulater::DelTriangulater():m_subdiv(MNull), m_triNum(0)
{
}

DelTriangulater::~DelTriangulater()
{
	uninit();
}

MRESULT DelTriangulater::doinit(MLong lWidth, MLong lHeight)
{
    MRECTF rectf = { -(MFloat)lWidth, -(MFloat)lHeight, (MFloat)lWidth*2, (MFloat)lHeight*2 };
//    MRECTF rectf = { (MFloat)0, (MFloat)0, (MFloat)lWidth, (MFloat)lHeight };
	m_subdiv = new Subdiv2D();
	m_subdiv->initDelaunay(rectf);

	return 0;
}

MRESULT DelTriangulater::uninit()
{
    if (NULL != m_subdiv) {
        delete(m_subdiv);
        m_subdiv = NULL;
    }
	return 0;
}

MRESULT DelTriangulater::setTriangulePoints(std::vector<MPOINTF> &points)
{
	for (std::vector<MPOINTF>::iterator it = points.begin(); it != points.end(); it++)
	{
		m_subdiv->insert(*it);
	}
	return 0;
}

MRESULT DelTriangulater::triVertexToIndex(std::vector<MPOINTF> &points, std::vector<MUInt32> &vertex)
{
	std::vector<VEC6F> triangle;
	m_subdiv->getTriangleList(triangle);
    m_triNum = (MInt32)triangle.size();
	for (int i = 0; i < triangle.size(); i++)
	{
		std::vector<MPOINTF> origTri;
		MPOINTF pt1 = { (MFloat)triangle[i].pt1x, (MFloat)triangle[i].pt1y };
		MPOINTF pt2 = { (MFloat)triangle[i].pt2x, (MFloat)triangle[i].pt2y };
		MPOINTF pt3 = { (MFloat)triangle[i].pt3x, (MFloat)triangle[i].pt3y };
		origTri.push_back(pt1);
		origTri.push_back(pt2);
		origTri.push_back(pt3);
		for (int j = 0; j < 3; j++)
		{
			for (size_t idx = 0; idx < points.size(); idx++)
			{
				if (points[idx].x == origTri[j].x && points[idx].y == origTri[j].y)
				{
					vertex.push_back((MInt32)idx);
					break;
				}
			}
		}
	}
	return 0;
}
