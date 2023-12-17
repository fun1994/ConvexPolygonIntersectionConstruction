#pragma once
#include <iostream>
#include "Segment.h"

class ConvexPolygonIntersectionConstruction {
	double area2(Point& p, Point& q, Point& r) {
		return p.x * q.y - p.y * q.x + q.x * r.y - q.y * r.x + r.x * p.y - r.y * p.x;
	}
	bool toLeft(Point& p, Point& q, Point& r) {
		return area2(p, q, r) > 0;
	}
	bool intersect(Segment& s1, Segment& s2) {
		return toLeft(*s1.first, *s1.second, *s2.first) != toLeft(*s1.first, *s1.second, *s2.second) && toLeft(*s2.first, *s2.second, *s1.first) != toLeft(*s2.first, *s2.second, *s1.second);
	}
	Point intersectionPoint(Segment& s1, Segment& s2) {
		double a1 = area2(*s1.first, *s2.first, *s2.second);
		double a2 = area2(*s1.second, *s2.first, *s2.second);
		Point p;
		p.x = a2 / (a2 - a1) * s1.first->x + a1 / (a1 - a2) * s1.second->x;
		p.y = a2 / (a2 - a1) * s1.first->y + a1 / (a1 - a2) * s1.second->y;
		return p;
	}
	bool inConvexPolygon(std::vector<Point>& P, Point& p) {
		for (int i = 0; i < P.size(); i++) {
			if (!toLeft(P[i], P[i == P.size() - 1 ? 0 : i + 1], p)) {
				return false;
			}
		}
		return true;
	}
	int LTL(std::vector<Point>& P) {
		int ltl = 0;
		for (int k = 1; k < P.size(); k++) {
			if (P[k].y < P[ltl].y || (P[k].y == P[ltl].y && P[k].x < P[ltl].x)) {
				ltl = k;
			}
		}
		return ltl;
	}
	int polarAnglePartition(std::vector<Point>& P, int low, int high) {
		Point pivot = P[low];
		while (low < high) {
			while (low < high && (toLeft(P[0], pivot, P[high]) || (P[high].x == pivot.x && P[high].y == pivot.y))) {
				high--;
			}
			P[low] = P[high];
			while (low < high && (toLeft(P[0], P[low], pivot) || (P[low].x == pivot.x && P[low].y == pivot.y))) {
				low++;
			}
			P[high] = P[low];
		}
		P[low] = pivot;
		return low;
	}
	void polarAngleQuickSort(std::vector<Point>& P, int low, int high) {
		if (low < high) {
			int pivot = polarAnglePartition(P, low, high);
			polarAngleQuickSort(P, low, pivot - 1);
			polarAngleQuickSort(P, pivot + 1, high);
		}
	}
	void polarAngleQuickSort(std::vector<Point>& P) {
		polarAngleQuickSort(P, 1, P.size() - 1);
	}
	double crossProduct(Segment& s1, Segment& s2) {
		return (s1.second->x - s1.first->x) * (s2.second->y - s2.first->y) - (s1.second->y - s1.first->y) * (s2.second->x - s2.first->x);
	}
	void advance(int& count, int& i, int n, bool flag, std::vector<Point>& P, Point& p) {
		count++;
		i = (i + 1) % n;
		if (flag) {
			P.push_back(p);
		}
	}
	int leftmost(std::vector<Point>& P) {
		int left = 0;
		int right = P.size() - 1;
		if (P[0].x < P[P.size() - 1].x) {
			while (left < right) {
				int mid = (left + right) / 2;
				if (P[mid].x < P[mid + 1].x) {
					right = mid;
				}
				else {
					if (P[mid].x > P[0].x) {
						right = mid - 1;
					}
					else {
						left = mid + 1;
					}
				}
			}
		}
		else {
			while (left < right) {
				int mid = (left + right) / 2;
				if (P[mid].x > P[mid + 1].x) {
					left = mid + 1;
				}
				else {
					if (P[mid].x > P[P.size() - 1].x) {
						left = mid + 1;
					}
					else {
						right = mid;
					}
				}
			}
		}
		return left;
	}
	int rightmost(std::vector<Point>& P) {
		int left = 0;
		int right = P.size() - 1;
		if (P[0].x < P[P.size() - 1].x) {
			while (left < right) {
				int mid = (left + right) / 2;
				if (P[mid].x < P[mid + 1].x) {
					left = mid + 1;
				}
				else {
					if (P[mid].x < P[P.size() - 1].x) {
						left = mid + 1;
					}
					else {
						right = mid;
					}
				}
			}
		}
		else {
			while (left < right) {
				int mid = (left + right) / 2;
				if (P[mid].x > P[mid + 1].x) {
					right = mid;
				}
				else {
					if (P[mid].x < P[0].x) {
						right = mid - 1;
					}
					else {
						left = mid + 1;
					}
				}
			}
		}
		return left;
	}
	void partition(std::vector<Point>& P, std::vector<Point>& PL, std::vector<Point>& PU, std::string index) {
		int left = leftmost(P);
		int right = rightmost(P);
		if (left < right) {
			PL.insert(PL.end(), P.begin() + left, P.begin() + (right + 1));
			PU.insert(PU.end(), P.rend() - (left + 1), P.rend());
			PU.insert(PU.end(), P.rbegin(), P.rend() - right);
		}
		else {
			PL.insert(PL.end(), P.begin() + left, P.end());
			PL.insert(PL.end(), P.begin(), P.begin() + (right + 1));
			PU.insert(PU.end(), P.rend() - (left + 1), P.rend() - right);
		}
		for (int i = 0; i < PL.size(); i++) {
			PL[i].type.push_back(index + "L");
			if (i == 0 || i == PL.size() - 1) {
				PL[i].type.push_back(index + "U");
			}
		}
		for (int i = 0; i < PU.size(); i++) {
			if (i == 0 || i == PU.size() - 1) {
				PU[i].type.push_back(index + "L");
			}
			PU[i].type.push_back(index + "U");
		}
	}
public:
	std::vector<Point> bruteForce(std::vector<Point>& P1, std::vector<Point>& P2) {
		std::vector<Point> P;
		for (int i = 0; i < P1.size(); i++) {
			if (inConvexPolygon(P2, P1[i])) {
				P.push_back(P1[i]);
			}
		}
		for (int j = 0; j < P2.size(); j++) {
			if (inConvexPolygon(P1, P2[j])) {
				P.push_back(P2[j]);
			}
		}
		for (int i = 0; i < P1.size(); i++) {
			Segment s1;
			s1.first = &P1[i];
			s1.second = &P1[i == P1.size() - 1 ? 0 : i + 1];
			for (int j = 0; j < P2.size(); j++) {
				Segment s2;
				s2.first = &P2[j];
				s2.second = &P2[j == P2.size() - 1 ? 0 : j + 1];
				if (intersect(s1, s2)) {
					Point p = intersectionPoint(s1, s2);
					P.push_back(p);
				}
			}
		}
		if (!P.empty()) {
			int ltl = LTL(P);
			std::swap(P[0], P[ltl]);
			polarAngleQuickSort(P);
		}
		return P;
	}
	std::vector<Point> ORourke(std::vector<Point>& P1, std::vector<Point>& P2) {
		std::vector<Point> P;
		int count1 = 0;
		int count2 = 0;
		int i1 = 0;
		int i2 = 0;
		int flag = 0;
		do {
			Segment s1, s2;
			s1.first = &P1[i1 % P1.size()];
			s1.second = &P1[(i1 + 1) % P1.size()];
			s2.first = &P2[i2 % P2.size()];
			s2.second = &P2[(i2 + 1) % P2.size()];
			bool flag1 = toLeft(*s2.first, *s2.second, *s1.second);
			bool flag2 = toLeft(*s1.first, *s1.second, *s2.second);
			if (intersect(s1, s2)) {
				if (flag == 0) {
					count1 = 0;
					count2 = 0;
				}
				Point p = intersectionPoint(s1, s2);
				P.push_back(p);
				if (flag1) {
					flag = 1;
				}
				else if (flag2) {
					flag = -1;
				}
			}
			if (crossProduct(s1, s2) > 0) {
				if (flag2) {
					advance(count1, i1, P1.size(), flag == 1, P, *s1.second);
				}
				else {
					advance(count2, i2, P2.size(), flag == -1, P, *s2.second);
				}
			}
			else {
				if (flag1) {
					advance(count2, i2, P2.size(), flag == -1, P, *s2.second);
				}
				else {
					advance(count1, i1, P1.size(), flag == 1, P, *s1.second);
				}
			}
		} while ((count1 < P1.size() || count2 < P2.size()) && count1 < 2 * P1.size() && count2 < 2 * P2.size());
		if (P.empty()) {
			if (inConvexPolygon(P2, P1[0])) {
				P = P1;
			}
			else if (inConvexPolygon(P1, P2[0])) {
				P = P2;
			}
		}
		return P;
	}
	std::vector<Point> planeSweeping(std::vector<Point>& P1, std::vector<Point>& P2) {
		std::vector<Point> P, PL, PU, P1L, P1U, P2L, P2U;
		partition(P1, P1L, P1U, "1");
		partition(P2, P2L, P2U, "2");
		int lower1 = -1;
		int upper1 = -1;
		int lower2 = -1;
		int upper2 = -1;
		Segment* s1L = NULL;
		Segment* s1U = NULL;
		Segment* s2L = NULL;
		Segment* s2U = NULL;
		Point* event = NULL;
		Point q;
		while (true) {
			Point* p = NULL;
			if ((event && event->x < P1L[0].x) || !event) {
				if ((p && P1L[0].x < p->x) || !p) {
					p = &P1L[0];
				}
			}
			if ((event && event->x < P2L[0].x) || !event) {
				if ((p && P2L[0].x < p->x) || !p) {
					p = &P2L[0];
				}
			}
			if (s1L && event->x < s1L->second->x) {
				if ((p && s1L->second->x < p->x) || !p) {
					p = s1L->second;
				}
			}
			if (s1U && event->x < s1U->second->x) {
				if ((p && s1U->second->x < p->x) || !p) {
					p = s1U->second;
				}
			}
			if (s2L && event->x < s2L->second->x) {
				if ((p && s2L->second->x < p->x) || !p) {
					p = s2L->second;
				}
			}
			if (s2U && event->x < s2U->second->x) {
				if ((p && s2U->second->x < p->x) || !p) {
					p = s2U->second;
				}
			}
			bool flag = false;
			Point r, s;
			if (s1L && s2L && intersect(*s1L, *s2L)) {
				r = intersectionPoint(*s1L, *s2L);
				r.type.push_back("1L");
				r.type.push_back("2L");
				if (event->x < r.x) {
					if ((p && r.x < p->x) || !p) {
						flag = true;
						s = r;
						p = &s;
					}
				}
			}
			if (s1L && s2U && intersect(*s1L, *s2U)) {
				r = intersectionPoint(*s1L, *s2U);
				r.type.push_back("1L");
				r.type.push_back("2U");
				if (event->x < r.x) {
					if ((p && r.x < p->x) || !p) {
						flag = true;
						s = r;
						p = &s;
					}
				}
			}
			if (s1U && s2L && intersect(*s1U, *s2L)) {
				r = intersectionPoint(*s1U, *s2L);
				r.type.push_back("1U");
				r.type.push_back("2L");
				if (event->x < r.x) {
					if ((p && r.x < p->x) || !p) {
						flag = true;
						s = r;
						p = &s;
					}
				}
			}
			if (s1U && s2U && intersect(*s1U, *s2U)) {
				r = intersectionPoint(*s1U, *s2U);
				r.type.push_back("1U");
				r.type.push_back("2U");
				if (event->x < r.x) {
					if ((p && r.x < p->x) || !p) {
						flag = true;
						s = r;
						p = &s;
					}
				}
			}
			if (flag) {
				q = s;
				p = &q;
			}
			if (p) {
				event = p;
			}
			else {
				break;
			}
			if (event->type.size() == 1) {
				if (event->type[0] == "1L") {
					if (s2L && s2U && toLeft(*s2L->first, *s2L->second, *event) && toLeft(*s2U->second, *s2U->first, *event)) {
						PL.push_back(*event);
					}
					lower1++;
					s1L->first = &P1L[lower1];
					s1L->second = &P1L[lower1 + 1];
				}
				else if (event->type[0] == "1U") {
					if (s2L && s2U && toLeft(*s2L->first, *s2L->second, *event) && toLeft(*s2U->second, *s2U->first, *event)) {
						PU.push_back(*event);
					}
					upper1++;
					s1U->first = &P1U[upper1];
					s1U->second = &P1U[upper1 + 1];
				}
				else if (event->type[0] == "2L") {
					if (s1L && s1U && toLeft(*s1L->first, *s1L->second, *event) && toLeft(*s1U->second, *s1U->first, *event)) {
						PL.push_back(*event);
					}
					lower2++;
					s2L->first = &P2L[lower2];
					s2L->second = &P2L[lower2 + 1];
				}
				else if (event->type[0] == "2U") {
					if (s1L && s1U && toLeft(*s1L->first, *s1L->second, *event) && toLeft(*s1U->second, *s1U->first, *event)) {
						PU.push_back(*event);
					}
					upper2++;
					s2U->first = &P2U[upper2];
					s2U->second = &P2U[upper2 + 1];
				}
			}
			else if (event->type.size() == 2) {
				if (event->type[0] == "1L" && event->type[1] == "1U") {
					if (s2L && s2U && toLeft(*s2L->first, *s2L->second, *event) && toLeft(*s2U->second, *s2U->first, *event)) {
						PL.push_back(*event);
						PU.push_back(*event);
					}
					if (event->x == P1L[0].x && event->y == P1L[0].y) {
						lower1++;
						s1L = new Segment;
						s1L->first = &P1L[lower1];
						s1L->second = &P1L[lower1 + 1];
						upper1++;
						s1U = new Segment;
						s1U->first = &P1U[upper1];
						s1U->second = &P1U[upper1 + 1];
					}
					else if (event->x == P1L[P1L.size() - 1].x && event->y == P1L[P1L.size() - 1].x) {
						lower1++;
						s1L = NULL;
						upper1++;
						s1U = NULL;
					}
				}
				else if (event->type[0] == "2L" && event->type[1] == "2U") {
					if (s1L && s1U && toLeft(*s1L->first, *s1L->second, *event) && toLeft(*s1U->second, *s1U->first, *event)) {
						PL.push_back(*event);
						PU.push_back(*event);
					}
					if (event->x == P2L[0].x && event->y == P2L[0].y) {
						lower2++;
						s2L = new Segment;
						s2L->first = &P2L[lower2];
						s2L->second = &P2L[lower2 + 1];
						upper2++;
						s2U = new Segment;
						s2U->first = &P2U[upper2];
						s2U->second = &P2U[upper2 + 1];
					}
					else if (event->x == P2L[P2L.size() - 1].x && event->y == P2L[P2L.size() - 1].x) {
						lower2++;
						s2L = NULL;
						upper2++;
						s2U = NULL;
					}
				}
				else if (event->type[0] == "1L" && event->type[1] == "2L") {
					PL.push_back(*event);
				}
				else if (event->type[0] == "1L" && event->type[1] == "2U") {
					PL.push_back(*event);
					PU.push_back(*event);
				}
				else if (event->type[0] == "1U" && event->type[1] == "2L") {
					PL.push_back(*event);
					PU.push_back(*event);
				}
				else if (event->type[0] == "1U" && event->type[1] == "2U") {
					PU.push_back(*event);
				}
			}
		}
		P.insert(P.end(), PL.begin(), PL.end());
		if (!PU.empty()) {
			P.insert(P.end(), PU.rbegin() + 1, PU.rend() - 1);
		}
		return P;
	}
};
