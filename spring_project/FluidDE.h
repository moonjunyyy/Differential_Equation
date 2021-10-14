#pragma once
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "bmpNew.h"

struct FluidMesh
{
	double rho = 997.; // kg/m3
	double Prs = 101325.; // Pa = kg / m s2
	double vscs = 1.79; // kg s / m2
	Eigen::Vector2d velo = { 0, 0 }, velo_p = { 0, 0 };
};
class FluidDE
{
public:
	FluidMesh *Mesh;
	int Sizex, Sizey;
	double dx = 0.01, dy = 0.01;
	double dt = 0.01;

	FluidDE(int x, int y)
	{
		Sizex = x, Sizey = y;
		Mesh = new FluidMesh[(Sizex + 1) * (Sizey + 1)];
	}
	~FluidDE()
	{
		delete[] Mesh;
	}

	double update()
	{
		FluidMesh* Cmesh;
		Cmesh = new FluidMesh[(Sizex + 1) * (Sizey + 1)];
		double Max = 0;

		for (int x = 0; x <= Sizex; x++)
		{
			for (int y = 0; y <= Sizey; y++)
			{
				Cmesh[x + y * (Sizex + 1)] = Mesh[x + y * (Sizex + 1)];
			}
		}

		for (int x = 1; x <= 500; x++)
		{
			for (int y = 1; y <= 500; y++)
			{
				double Xp	= Cmesh[((x - 1) < 1 ? 0 : (x - 1)) + (y * (Sizex + 1))].velo_p(0),
					Xpp		= Cmesh[((x - 2) < 1 ? 0 : (x - 2)) + (y * (Sizex + 1))].velo_p(0),
					Yp		= Cmesh[ x + ((y - 1) < 1 ? 0 : (y - 1)) * (Sizex + 1) ].velo_p(1),
					Ypp		= Cmesh[ x + ((y - 2) < 1 ? 0 : (y - 2)) * (Sizex + 1) ].velo_p(1);

				FluidMesh& This = Mesh[x + y * (Sizex + 1)];
				This.velo = ((This.vscs / This.rho) * Eigen::Vector2d( (This.velo_p(0) + Xpp - Xp - Xp) / (dx * dx), (This.velo_p(1) + Ypp - Yp - Yp) / (dy * dy))
					- Eigen::Vector2d( This.velo_p(0) * (This.velo_p(0) - Xp) / dx, This.velo_p(1) * (This.velo_p(1) - Yp) / dy)) * dt + This.velo_p;

				This.velo_p = This.velo;
				if (This.velo_p.norm() > Max) { Max = This.velo_p.norm(); }
			}
		}
		delete[] Cmesh;
		return Max;
	}
	double DrawIMG(char* FileName)
	{
		double Max = update();
		unsigned char* IMG = new unsigned char[Sizex * Sizey * 3];
		for (int x = 0; x < 500; x++)
		{
			for (int y = 0; y < 500; y++)
			{
				double abv = Mesh[(x + 1) + (y + 1) * (Sizex + 1)].velo_p.norm() * 100000 * 255;
				IMG[(x + y * Sizex) * 3 + 0] = abv > 255 ? 255 : (unsigned char)abv;
				IMG[(x + y * Sizex) * 3 + 1] = abv > 255 ? 255 : (unsigned char)abv;
				IMG[(x + y * Sizex) * 3 + 2] = abv > 255 ? 255 : (unsigned char)abv;
			}
		}
		WriteBmp(FileName, IMG, Sizex, Sizey);
		delete[] IMG;
		return Max;
	}
};