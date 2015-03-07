/*  DeClone: A software for computing and analyzing ancestral adjacency scenarios.
 *  Copyright (C) 2015 Cedric Chauve, Yann Ponty, Ashok Rajaraman, Joao P.P. Zanetti
 *
 *  This file is part of DeClone.
 *  
 *  DeClone is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DeClone is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DeClone.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Contact: <yann.ponty@lix.polytechnique.fr>.
 *
 *
 *  DeClone uses the Quickhull algorithm implementation programmed by 
 *  Anatoly V. Tomilov. The code is available on <https://bitbucket.org/tomilov/quickhull/src/585267abb3a63794c04fc8325aa9ec9f726112ed/include/quickhull.hpp?at=master>.
 *
 *  Contact: <tomilovanatoliy@gmail.com>
 */

#ifndef SVGDRIVER_HH
#define SVGDRIVER_HH

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Point{
  public:
		double x,y;
		Point(const Point & p){
			x=p.x;
			y=p.y;
		}
		Point(double a, double b){
			x=a;
			y=b;
		}
};

Point operator+ (const Point & p1, const Point & p2);
Point operator* (const Point & p,double k);
Point operator* (double k, const Point & p);


class Color{
	public:
		double r,g,b,a;
		
		Color(){
			r = 0.;
			g = 0.;
			b = 0.;
			a = 1.;
		};


		Color(const Color & c){
			r = c.r;
			g = c.g;
			b = c.b;
			a = 1.;
		};

		Color(double ar, double ag, double ab){
			r = ar;
			g = ag;
			b = ab;
			a = 1.;
		};

		Color(double ar, double ag, double ab, double aa){
			r = ar;
			g = ag;
			b = ab;
			a = aa;
		};
		
		double getAlpha()
		{return a;}

};

ostream & operator<< (ostream & o, const Color & c);

class SVGFile{
	private:
		ofstream out;
		Color color;
		double thickness;
	public:
		
		SVGFile(string path) : out(path.c_str()){
			thickness = 1.0;
			color = Color(0.,0.,0.);
			out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
				<< "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n"
				<< "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
				<< "\n"
				<< "<svg width=\"100%\" height=\"100%\" version=\"1.1\"\n"
				<< "xmlns=\"http://www.w3.org/2000/svg\">\n";
		}
		
		~SVGFile(){
			out << "</svg>\n";
			out.close();
		}

		void setColor(const Color & c)
		{
			color = c;
		}
		
		void drawLine(Point p1, Point p2)
		{
				out << "<line x1=\"" << p1.x << "\" y1=\"" << p1.y
				<< "\" x2=\"" << p2.x << "\" y2=\"" << p2.y
				<< "\" stroke=\"" << color
				<< "\" stroke-width=\"" << thickness << "opacity:"<< color.getAlpha()<<";\" />\n";
		}

		void drawString(Point p, string s, double fontsize)
		{
			double x = p.x;
			double y = (p.y+(fontsize/3.));
			
				out << "<text x=\""
				<< x
				<< "\" y=\""
				<< y
				<< "\" text-anchor=\"middle\" font-family=\"Verdana\" font-size=\""
				<< fontsize << "\" fill=\"" << color << "\" >"
				<< s << "</text>\n";
		}

		void drawStringLeftAlign(Point p, string s, double fontsize)
		{
			double x = p.x;
			double y = (p.y+(fontsize/3.));
			
				out << "<text x=\""
				<< x
				<< "\" y=\""
				<< y
				<< "\" text-anchor=\"left\" font-family=\"Verdana\" font-size=\""
				<< fontsize << "\" fill=\"" << color << "\" >"
				<< s << "</text>\n";
		}

		void drawStringRotated90(Point p, string s, double fontsize)
		{
			double x = p.x-(fontsize/3.);
			double y = (p.y);
			
				out << "<text x=\""
				<< x
				<< "\" y=\""
				<< y
				<< "\" text-anchor=\"middle\" font-family=\"Verdana\" font-size=\""
				<< fontsize << "\" fill=\"" << color << "\" transform=\"rotate(90 "<<x<<" "<<y<<")\" >"
				<< s << "</text>\n";
		}

		void fillSquare(Point center, double side)
		{
			vector<Point> ps;
			ps.push_back(Point(center.x-side/2.,center.y-side/2.));
			ps.push_back(Point(center.x+side/2.,center.y-side/2.));
			ps.push_back(Point(center.x+side/2.,center.y+side/2.));
			ps.push_back(Point(center.x-side/2.,center.y+side/2.));
			ps.push_back(Point(center.x-side/2.,center.y-side/2.));
			fillPolygon(ps);
		}
		
		void drawSquare(Point center, double side)
		{
			vector<Point> ps;
			ps.push_back(Point(center.x-side/2.,center.y-side/2.));
			ps.push_back(Point(center.x+side/2.,center.y-side/2.));
			ps.push_back(Point(center.x+side/2.,center.y+side/2.));
			ps.push_back(Point(center.x-side/2.,center.y+side/2.));
			ps.push_back(Point(center.x-side/2.,center.y-side/2.));
			drawPolygon(ps);
		}

		void drawCircle(Point center, double radius)
		{
			out << "<circle cx=\""<< center.x<<"\" cy=\""<< center.y<<"\" r=\""<< radius
			<< "\" style=\"fill:none; stroke:" << color << "; stroke-width:" << thickness << "; opacity:"<< color.getAlpha()<<";\" />\n";
		}

		void fillCircle(Point center, double radius)
		{
			out << "<circle cx=\""<< center.x<<"\" cy=\""<< center.y<<"\" r=\""<< radius
			<< "\" style=\"fill:" << color << "; stroke:none; stroke-width:" << thickness << ";opacity:"<< color.getAlpha()<<";\" />\n";
		}
		

		void drawPolygon(vector<Point> points) {
			out << "<path d=\"";
			for (int i = 0; i < points.size(); i++) {
				if (i == 0) {
					out << "M " << points[i].x << " " << points[i].y
							<< " ";
				} else {
					out << "L " << points[i].x << " " << points[i].y
							<< " ";
				}
			}
			out << "z\" style=\"fill:none; stroke:" << color << "; stroke-width:" << thickness << ";opacity:"<< color.getAlpha()<<";\"/>\n";
		}

		void fillPolygon(vector<Point> points) {
			out << "<path d=\"";
			for (int i = 0; i < points.size(); i++) {
				if (i == 0) {
					out << "M " << points[i].x << " " << points[i].y
							<< " ";
				} else {
					out << "L " << points[i].x << " " << points[i].y
							<< " ";
				}
			}
			out << "z\" style=\"fill:" << color << ";stroke:none; stroke-width:" << thickness << ";opacity:"<< color.getAlpha()<<";\"/>\n";
		}
};

#endif
