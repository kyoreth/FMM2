// FMM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include "octree.h"

std::ifstream infile("charges.txt");

//globals

int leveladd = 0;
int maxlvl = 0;
double probe_pot;
double probe_force_x = 0, probe_force_y = 0, probe_force_z = 0;


struct point_charge {
	double q, x, y, z;
};

struct moment_list {
	double Q0, QX, QY, QZ, QXX, QXY, QXZ, QYY, QZZ, QYZ;
};

struct tcell {
	int lvl = 0;
	double x, y, z;
	double potential;
	double force_x, force_y, force_z;
};

std::vector<tcell> tree;

moment_list moment(std::vector<point_charge>& points, double cx, double cy, double cz) {

	int N = points.size();

	double Q0 = 0, QX = 0, QY = 0, QZ = 0, QXX = 0, QXY = 0, QXZ = 0, QYY = 0, QZZ = 0, QYZ = 0;


	for (int i = 0; i < N; i++) {

		double x = points[i].x - cx;
		double y = points[i].y - cy;
		double z = points[i].z - cz;
		double q = points[i].q;

		Q0 += q;
		QX += x*q;
		QY += y*q;
		QZ += z*q;
		QXX += pow(x,2)*q;
		QYY += pow(y,2)*q;
		QZZ += pow(z,2)*q;
		QXY += x*y*q;
		QXZ += x*z*q;
		QYZ += y*z*q;

	}

	QXX *= 1 / 2;
	QYY *= 1 / 2;
	QZZ *= 1 / 2;
	QXY *= 1 / 2;
	QXZ *= 1 / 2;
	QYZ *= 1 / 2;

	moment_list moments = moment_list{ Q0, QX, QY, QZ, QXX, QXY, QXZ, QYY, QZZ, QYZ };

	return(moments);

}


//given a set of point charges, create a tree of cells and calculate each cells moment and potential
void create_tree(std::vector<point_charge>& points) {

	int N = points.size();

	if (N == 0) {
		leveladd = -1;
		return;
	}

	tcell cell;

	cell.lvl += leveladd;

	if (cell.lvl > maxlvl) { maxlvl = cell.lvl; };

	leveladd = 1;
	
	double sum_x = 0, sum_y = 0, sum_z = 0;

	for (int i = 0; i < N; i++) {
		sum_x += points[i].x;
		sum_y += points[i].y;
		sum_z += points[i].z;
	}

	cell.x = sum_x / N;
	cell.y = sum_y / N;
	cell.z = sum_z / N;

	moment_list moments = moment(points, cell.x, cell.y, cell.z);

	double Q0 = moments.Q0;
	double QX = moments.QX;
	double QY = moments.QY;
	double QZ = moments.QZ;
	double QXX = moments.QXX;
	double QYY = moments.QYY;
	double QZZ = moments.QZZ;
	double QXY = moments.QXY;
	double QXZ = moments.QXZ;
	double QYZ = moments.QYZ;


	double U = cell.x / (pow(cell.x, 2) + pow(cell.y, 2) + pow(cell.z, 2));

	double potential_QXX = QXX*(((pow(U, (3 / 2)) - (3 * pow(cell.x, 2) * sqrt(U))) / pow(U, 3)));
	double potential_QYY = QYY*(((pow(U, (3 / 2)) - (3 * pow(cell.y, 2) * sqrt(U))) / pow(U, 3)));
	double potential_QZZ = QZZ*(((pow(U, (3 / 2)) - (3 * pow(cell.z, 2) * sqrt(U))) / pow(U, 3)));
	double potential_QXY = QXY*-(3 * cell.x * cell.y / pow(U, (5 / 2)));
	double potential_QXZ = QXZ*-(3 * cell.x * cell.z / pow(U, (5 / 2)));
	double potential_QYZ = QYZ*-(3 * cell.y * cell.z / pow(U, (5 / 2)));

	//electric field components

	//n=000
	double force_x_0 = Q0 * (cell.x / pow(U, (3 / 2)));
	double force_y_0 = Q0 * (cell.y / pow(U, (3 / 2)));
	double force_z_0 = Q0 * (cell.z / pow(U, (3 / 2)));

	//n=100
	double force_x_100 = -QX * ((pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2))) / pow(U, 3));
	double force_y_100 = -QX * ((3 * cell.y*cell.x) / pow(U, (5 / 2)));
	double force_z_100 = -QX * ((3 * cell.z*cell.x) / pow(U, (5 / 2)));

	//n=010
	double force_y_010 = -QY * ((pow(U, (3 / 2)) - 3 * pow(cell.y, 2)*pow(U, (1 / 2))) / pow(U, 3));
	double force_x_010 = -QY * ((3 * cell.y*cell.x) / pow(U, (5 / 2)));
	double force_z_010 = -QY * ((3 * cell.z*cell.y) / pow(U, (5 / 2)));

	//n=001
	double force_z_001 = -QZ * ((pow(U, (3 / 2)) - 3 * pow(cell.z, 2)*pow(U, (1 / 2))) / pow(U, 3));
	double force_x_001 = -QZ * ((3 * cell.z*cell.x) / pow(U, (5 / 2)));
	double force_y_001 = -QZ * ((3 * cell.z*cell.y) / pow(U, (5 / 2)));

	//n=200
	double force_x_200 = -QXX *(((-6 * pow(cell.x, 3) - 3 * cell.x*pow(cell.z, 2) - 3 * cell.x*pow(cell.y, 2))*pow(U, (5 / 2)) - 6 * cell.x*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	double force_y_200 = -QXX *(((3 * pow(cell.y, 3) + 3 * cell.y*pow(cell.z, 2))*pow(U, (5 / 2)) - 6 * cell.y*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	double force_z_200 = -QXX *(((3 * pow(cell.z, 3) + 3 * cell.z*pow(cell.y, 2))*pow(U, (5 / 2)) - 6 * cell.z*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	
	//n=020
	double force_y_020 = -QYY *(((-6 * pow(cell.y, 3) - 3 * cell.y*pow(cell.z, 2) - 3 * cell.y*pow(cell.x, 2))*pow(U, (5 / 2)) - 6 * cell.y*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	double force_x_020 = -QYY *(((3 * pow(cell.x, 3) + 3 * cell.x*pow(cell.z, 2))*pow(U, (5 / 2)) - 6 * cell.x*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	double force_z_020 = -QYY *(((3 * pow(cell.z, 3) + 3 * cell.z*pow(cell.y, 2))*pow(U, (5 / 2)) - 6 * cell.z*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));

	//n=002
	double force_z_002 = -QZZ *(((-6 * pow(cell.z, 3) - 3 * cell.z*pow(cell.y, 2) - 3 * cell.z*pow(cell.x, 2))*pow(U, (5 / 2)) - 6 * cell.z*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	double force_x_002 = -QZZ *(((3 * pow(cell.x, 3) + 3 * cell.x*pow(cell.z, 2))*pow(U, (5 / 2)) - 6 * cell.x*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));
	double force_y_002 = -QZZ *(((3 * pow(cell.y, 3) + 3 * cell.y*pow(cell.z, 2))*pow(U, (5 / 2)) - 6 * cell.y*pow(U, 2)*(pow(U, (3 / 2)) - 3 * pow(cell.x, 2)*pow(U, (1 / 2)))) / pow(U, 6));

	//n=110
	double force_x_110 = -QXY * 3 * cell.y*((pow(U, (5 / 2)) - 5 * pow(cell.x, 2)*pow(U, (3 / 2))) / pow(U, 5));
	double force_y_110 = -QXY * 3 * cell.x*((pow(U, (5 / 2)) - 5 * pow(cell.y, 2)*pow(U, (3 / 2))) / pow(U, 5));
	double force_z_110 = -QXZ * 15 * ((cell.x*cell.y*cell.z) / pow(U, (7 / 2)));

	//n=101
	double force_x_101 = -QXY * 3 * cell.z*((pow(U, (5 / 2)) - 5 * pow(cell.x, 2)*pow(U, (3 / 2))) / pow(U, 5));
	double force_z_101 = -QXY * 3 * cell.x*((pow(U, (5 / 2)) - 5 * pow(cell.z, 2)*pow(U, (3 / 2))) / pow(U, 5));
	double force_y_101 = -QXZ * 15 * ((cell.x*cell.y*cell.z) / pow(U, (7 / 2)));

	//n=011
	double force_z_011 = -QYZ * 3 * cell.y*((pow(U, (5 / 2)) - 5 * pow(cell.z, 2)*pow(U, (3 / 2))) / pow(U, 5));
	double force_y_011 = -QYZ * 3 * cell.z*((pow(U, (5 / 2)) - 5 * pow(cell.y, 2)*pow(U, (3 / 2))) / pow(U, 5));
	double force_x_011 = -QYZ * 15 * ((cell.x*cell.y*cell.z) / pow(U, (7 / 2)));

	cell.force_x = force_x_0 + force_x_100 + force_x_010 + force_x_001 + force_x_200 + force_x_020 + force_x_002 + force_x_110 + force_x_101 + force_x_011;
	cell.force_y = force_y_0 + force_y_100 + force_y_010 + force_y_001 + force_y_200 + force_y_020 + force_y_002 + force_y_110 + force_y_101 + force_y_011;
	cell.force_x = force_z_0 + force_z_100 + force_z_010 + force_z_001 + force_z_200 + force_z_020 + force_z_002 + force_z_110 + force_z_101 + force_z_011;


	cell.potential = Q0 + QX*(cell.x / pow(U, (3 / 2))) + QY*(cell.y / pow(U, (3 / 2))) + QZ*(cell.z / pow(U, (3 / 2)));
	
	//+potential_QXX + potential_QYY + potential_QZZ + potential_QXY + potential_QXZ + potential_QYZ);
	


	tree.insert(tree.end(), cell);

	std::vector<point_charge> a, b, c, d, e, f, g, h;

	for (int i = 0; i < N; i++) {
		if (points[i].x < cell.x && points[i].y < cell.y && points[i].z < cell.z) {
			a.insert(a.end(), points[i]);
		}

		if (points[i].x > cell.x && points[i].y < cell.y && points[i].z < cell.z) {
			b.insert(b.end(), points[i]);
		}

		if (points[i].x < cell.x && points[i].y > cell.y && points[i].z < cell.z) {
			c.insert(c.end(), points[i]);
		}

		if (points[i].x > cell.x && points[i].y > cell.y && points[i].z < cell.z) {
			d.insert(d.end(), points[i]);
		}

		if (points[i].x < cell.x && points[i].y < cell.y && points[i].z > cell.z) {
			e.insert(e.end(), points[i]);
		}

		if (points[i].x > cell.x && points[i].y < cell.y && points[i].z > cell.z) {
			f.insert(f.end(), points[i]);
		}

		if (points[i].x < cell.x && points[i].y > cell.y && points[i].z > cell.z) {
			g.insert(g.end(), points[i]);
		}

		if (points[i].x > cell.x && points[i].y > cell.y && points[i].z > cell.z) {
			h.insert(h.end(), points[i]);
		}
	}

	create_tree(a);
	create_tree(b);
	create_tree(c);
	create_tree(d);
	create_tree(e);
	create_tree(f);
	create_tree(g);
	create_tree(h);

}




int main()
{




	std::vector<point_charge> PC_list;

	double x, y, z, q;


	while (infile >> x >> y >> z >> q)
	{
		point_charge PC;
		PC.x = x;
		PC.y = y;
		PC.z = z;
		PC.q = q;

		PC_list.insert(PC_list.end(), PC);
	}

	std::cout << PC_list.size() << std::endl;

	create_tree(PC_list);

	//probe location
	std::vector<int> probe{ 0,0,0 };

	//distance considered far = outer limit
	int F = 300;
	
	//how far the far distance creeps closer each step
	int J = 50;

	int N = tree.size();

	//calculate potential att origo 
	
	/*for (int i = maxlvl; i >= 0; i--) {
		for (int n = 0; n < N; n++) {
			if ((tree[n].lvl == i) && (F > sqrt(pow(tree[n].x, 2) + pow(tree[n].y, 2) + pow(tree[n].y, 2)))) {
				probe_pot += tree[n].potential;
				probe_force_x += tree[n].force_x;
				probe_force_y += tree[n].force_y;
				probe_force_z += tree[n].force_z;
			}
		}

		F += J;
	}*/

	probe_pot = tree[0].potential;
	probe_force_x = tree[0].force_x;
	probe_force_y = tree[0].force_y;
	probe_force_z = tree[0].force_z;

	std::cout << probe_pot << std::endl << "Force X:" << probe_force_x <<" Y:" << probe_force_y << " Z:" << probe_force_z << std::endl;

	return 0;
}

