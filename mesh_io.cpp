#include "stdafx.h"

#include "mesh.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

struct F2CIndex
{
	std::istream *s;

	friend F2CIndex& operator>>(std::istream& is, F2CIndex &i);

	operator std::istream&() {
		return *s;
	}

	F2CIndex& operator>>(int &i) {
		*s >> i;
		i -= 1;
		return *this;
	}
};
F2CIndex& operator>>(std::istream& is, F2CIndex &i) {
	i.s = &is;
	return i;
}


enum VtkCellType
{
	VTK_LINE = 3,
	VTK_TRIANGLE = 5,
	VTK_QUAD = 9,
	VTK_TETRA = 10,
	VTK_HEXAHEDRON = 12,
	VTK_QUADRATIC_EDGE = 21,
	VTK_QUADRATIC_TRIANGLE = 22,
	VTK_QUADRATIC_QUAD = 23,
	VTK_QUADRATIC_TETRA = 24,

};


template<typename T>
void writeBigEndian(FILE* fp, const T* data, size_t counts)
{
	const char* p = (const char*)data;
	char buffer[sizeof(T) * 256];
	for (size_t written = 0; written<counts; written += _countof(buffer) / sizeof(T)) {
		size_t toWrite = std::min(_countof(buffer) / sizeof(T), counts - written);

		for (int i = 0; i<toWrite; i++) {
			for (int j = 0; j<sizeof(T); j++) {
				//swap bytes
				buffer[i * sizeof(T) + j] = p[sizeof(T)*(i + 1) - j - 1];
			}
		}

		fwrite(buffer, sizeof(T) * toWrite, 1, fp);
	}
}

inline void writeBigEndianCastFloat(FILE* fp, const double* data, size_t count)
{
	float buffer[256];
	size_t index = 0;
	for (;;) {
		size_t packed = 0;
		for (size_t i = 0; i<_countof(buffer) && i + index<count; i++, packed++) {
			buffer[i] = (float)data[index + i];
		}
		writeBigEndian(fp, buffer, packed);
		if (packed<_countof(buffer))break;
		index += packed;
	}
}

////////////////////////////////////////////////////////////////////////////////



void Mesh::printInfo(std::ostream & os)
{
	os << "Mesh:" << std::endl;
	os << "nvert=" << numVertices() << std::endl;
	os << "nnode=" << numNodes() << std::endl;
	os << "nelem=" << numElements() << std::endl;
	os << "nbndry=" << numBoundaryElements() << std::endl;
}

void Mesh::build2DSquare(const int ncell[], const double len[], const double domlo[])
{
	std::cout << __FUNCTION__ << std::endl;

	//
	const int sdim = sfem::SpaceDim();
	assert(sdim == 2);

	//
	const int nx = ncell[0];
	const int ny = ncell[1];
	const double lx = len[0];
	const double ly = len[1];
	const double xlo = domlo[0];
	const double ylo = domlo[1];
	const double dx = lx / nx;
	const double dy = ly / ny;

	// index function
	auto index = [nx](int i, int j) -> int {
		return i + j * (nx + 1);
	};

	// nodes
	for (int j = 0; j <= ny; j++) {
		double y = ylo + dy * j;
		for (int i = 0; i <= nx; i++) {
			double x = xlo + dx * i;
			VectorXd pos(sdim);
			pos << x, y;
			nodes.push_back(pos);
		}
	}

	// elements
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			int ind0 = index(i, j);
			int ind1 = index(i + 1, j);
			int ind2 = index(i + 1, j + 1);
			int ind3 = index(i, j + 1);

			Element *elem = nullptr;
			if (1) { // two triangles each square
				elem = new TriangleElement(ind0, ind1, ind2);
				elem->setAttribute(0);
				elements.push_back(elem);

				elem = new TriangleElement(ind0, ind2, ind3);
				elem->setAttribute(0);
				elements.push_back(elem);
			}
			else { // one quadrilateral each square
				//elem = new 
				// TODO
			}
		}
	}

	// boundary segments
	for (int i = 0; i < nx; i++) {
		int ind0 = index(i, 0);
		int ind1 = index(i + 1, 0);
		Element *e = new SegmentElement(ind0, ind1);
		e->setAttribute(1);
		boundary.push_back(e);
	}
	for (int j = 0; j < ny; j++) {
		int ind0 = index(nx, j);
		int ind1 = index(nx, j + 1);
		Element *e = new SegmentElement(ind0, ind1);
		e->setAttribute(2);
		boundary.push_back(e);
	}
	for (int i = nx - 1; i >= 0; i--) {
		int ind0 = index(i + 1, ny);
		int ind1 = index(i, ny);
		Element *e = new SegmentElement(ind0, ind1);
		e->setAttribute(3);
		boundary.push_back(e);
	}
	for (int j = ny - 1; j >= 0; j--) {
		int ind0 = index(0, j + 1);
		int ind1 = index(0, j);
		Element *e = new SegmentElement(ind0, ind1);
		e->setAttribute(4);
		boundary.push_back(e);
	}

}

////////////////////////////////////////////////////////////////////////////////


int Mesh::readGmsh(const char * filename)
{
	std::cout << __FUNCTION__ << ": " << filename << std::endl;

	// space dimensions
	const int sdim = ::sfem::SpaceDim();
	assert(sdim == 2 || sdim == 3);

	// file stream
	std::ifstream fs(filename);

	// line buffer
	std::string line;

	// physical names
	std::vector<std::pair<int, std::string> > phys_names;

	// Section: MeshFormat
	{
		const std::string sec = "MeshFormat";
		std::getline(fs, line); assert(line == "$" + sec);
		std::getline(fs, line);
		std::getline(fs, line); assert(line == "$End" + sec);
	}
	
	// Section: PhysicalNames
	{
		const std::string sec = "PhysicalNames";
		std::getline(fs, line); assert(line == "$" + sec);

		//
		std::getline(fs, line);
		const int num = std::stoi(line);
		std::cout << sec << ": n=" << num << std::endl;

		// 
		for (int i = 0; i < num; i++) {
			std::getline(fs, line);
			std::istringstream iss(line);

			int nd; iss >> nd;
			//int idx; iss >> idx; idx -= 1;
			int idx; iss >> F2CIndex() >> idx;
			std::string name; iss >> name;
			// remove quotes on name
			if (!name.empty()) {
				name.erase(name.begin());
				name.erase(name.end()-1);
			}

			std::cout << "idx=" << idx << "; name=" << name << "; dim=" << nd << std::endl;
			
			phys_names.push_back(std::make_pair(nd, name));
		}

		std::getline(fs, line); assert(line == "$End" + sec);
	}

	// Section: Nodes
	{
		const std::string sec = "Nodes";
		std::getline(fs, line); assert(line == "$" + sec);

		//
		std::getline(fs, line);
		const int num = std::stoi(line);
		std::cout << sec << ": n=" << num << std::endl;

		for (int i = 0; i < num; i++) {
			std::getline(fs, line);
			std::istringstream iss(line);

			int idx;
			VectorXd vpos(sdim);
			iss >> idx;
			for (int dim = 0; dim < sdim; dim++) {
				iss >> vpos(dim);
			}

			// save node
			nodes.push_back(vpos);
		}

		std::getline(fs, line); assert(line == "$End" + sec);
	}

	// Section: Elements
	{
		const std::string sec = "Elements";
		std::getline(fs, line); assert(line == "$" + sec);

		//
		std::getline(fs, line);
		const int num = std::stoi(line);
		std::cout << sec << ": n=" << num << std::endl;

		for (int i = 0; i < num; i++) {
			std::getline(fs, line);
			std::istringstream iss(line);

			int idx; iss >> idx;
			int elem_type; iss >> elem_type;
			int num_tag; iss >> num_tag;
			int phys_idx; iss >> F2CIndex() >> phys_idx;
			for (int k = 1; k < num_tag; k++) {
				int tmp; iss >> tmp;
			}

			// 
			Element *elem = nullptr;

			if (elem_type == 1) { // line
				int nodes[2];
				for (int k = 0; k < _countof(nodes); k++) {
					iss >> F2CIndex() >> nodes[k];
				}
				elem = new SegmentElement(nodes);
			}
			else if (elem_type == 2) { // triangle
				int nodes[3];
				for (int k = 0; k < _countof(nodes); k++) {
					iss >> F2CIndex() >> nodes[k];
				}
				elem = new TriangleElement(nodes);
			}
			else if (elem_type == 8) { // 2nd-order line 
				int nodes[3];
				for (int k = 0; k < _countof(nodes); k++) {
					iss >> F2CIndex() >> nodes[k];
				}
				elem = new QuadraticSegmentElement(nodes);
			}
			else if (elem_type == 9) { // 2nd-order triangle
				int nodes[6];
				for (int k = 0; k < _countof(nodes); k++) {
					iss >> F2CIndex() >> nodes[k];
				}
				elem = new QuadraticTriangleElement(nodes);
			}
			else {
				std::cerr << "Unsupported element type=" << elem_type << std::endl;
				return 1;
			}

			elem->setAttribute(phys_idx);

			// the physical group
			const auto& phys = phys_names[phys_idx];
			if (phys.first < sdim) { 
				// low dimension, boundary element
				this->boundary.push_back(elem);
			}
			else { 
				// domain element
				this->elements.push_back(elem);
			}

		}

		std::getline(fs, line); assert(line == "$End" + sec);
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////


int Mesh::writeVtk(const char * filename)
{
	std::cout << __FUNCTION__ << ": " << filename << std::endl;

	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		std::cerr << "Failed to write " << filename << std::endl;
		return 1;
	}

	//
	const int sdim = ::sfem::SpaceDim();

	//header
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	//title
	fprintf(fp, "SFEM mesh\n");
	//type
	fprintf(fp, "BINARY\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	//
	// field data here
	//

	//
	// nodes
	//
	const int num_node = numNodes();
	
	fprintf(fp, "POINTS %d float\n", num_node);
	for (int i = 0; i < num_node; i++) {
		double data[3] = { 0 };
		std::copy(nodes[i].data(), nodes[i].data() + sdim, data);
		writeBigEndianCastFloat(fp, data, 3);
	}
	fprintf(fp, "\n");

	//
	// cells
	//
	const int num_elem = numElements();

	int count = num_elem;
	for (const Element *e : elements) {
		count += e->numNodes();
	}

	fprintf(fp, "CELLS %d %d\n", num_elem, count);
	for (int i = 0; i < num_elem; i++) {
		int nn = elements[i]->numNodes();
		std::vector<int> data(nn + 1);
		data[0] = nn;
		for (int k = 1; k <= nn; k++) {
			data[k] = elements[i]->getIndices()[k - 1];
		}
		writeBigEndian(fp, data.data(), nn + 1);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", num_elem);
	for (int i = 0; i < num_elem; i++) {
		int tt = elements[i]->type();
		int data = -1;
		if (tt == ElementType::Tri) {
			data = VTK_TRIANGLE;
		}
		else if (tt == ElementType::Quad) {
			data = VTK_QUAD;
		}
		else if (tt == ElementType::Tri2) {
			data = VTK_QUADRATIC_TRIANGLE;
		}
		else {
			std::cerr << "Unsupported element type=" << tt << std::endl;
			return 1;
		}

		writeBigEndian(fp, &data, 1);
	}
	fprintf(fp, "\n");


	fclose(fp);

	return 0;
}



////////////////////////////////////////////////////////////////////////////////


void VisDataSet::writeVtk()
{
	std::string filename = title + ".vtk";

	mesh->writeVtk(filename.c_str());

	// append data
	FILE *fp = fopen(filename.c_str(), "ab");

	fprintf(fp, "POINT_DATA %d\n", mesh->numNodes());

	for (auto it = dataset.cbegin(); it != dataset.cend(); ++it) {
		fprintf(fp, "SCALARS %s float 1\n", it->first.c_str());
		fprintf(fp, "LOOKUP_TABLE default\n");
		const auto &data = *it->second;
		writeBigEndianCastFloat(fp, data.data(), mesh->numNodes());
		fprintf(fp, "\n");
	}

	fclose(fp);
}



////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

