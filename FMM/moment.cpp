#include "stdafx.h"
#include <array>

void moment(points) {

	int N = points.size()/4;

	double Q0, QX, QY, QZ, QXX, QXY, QXZ, QYY, QZZ, QYZ = 0;

	for (i = 0; i <= N; i++) {

		Q0 += points[i][0];
		QX += points[i][1];
		QY += points[i][2];
		QZ += points[i][3];
		QXX += points[i][1] * points[i][1];
		QXY += points[i][2] * points[i][1];
		QXZ += points[i][3] * points[i][1];
		QYZ += points[i][2] * points[i][3];
		QYY += points[i][2] * points[i][2];
		QZZ += points[i][3] * points[i][3];

	}
}