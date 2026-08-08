#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <iostream>
#include <trackcpp/output.h>
#include <trackcpp/optics.h>
#include <trackcpp/dynap.h>
#include <trackcpp/tracking.h>
#include <trackcpp/lattice.h>
#include <trackcpp/flat_file.h>
#include <trackcpp/kicktable.h>
#include <trackcpp/passmethods.h>
#include <trackcpp/accelerator.h>
#include <trackcpp/elements.h>
#include <trackcpp/pos.h>
#include <trackcpp/tpsa.h>
#include <trackcpp/auxiliary.h>
#include <trackcpp/multithreads.h>
#include <trackcpp/linalg.h>
#include <iostream>

namespace py = pybind11;

py::array_t<double> get_double_ptr(const double *double_ptr, const unsigned int s1, const unsigned int s2 = 0)
{
    if (s2 == 0)
    {
        return py::array_t<double>(
            {s1},             // Shape of the array
            {sizeof(double)}, // Strides of the array (assuming it's a 1D array)
            double_ptr,       // Pointer to the data
            py::capsule(double_ptr, [](void *double_ptr)
                        {
                            // Deleter function, called when the Python object is destroyed
                            // Do nothing here because we don't own the memory
                        }));
    }
    else
    {
        return py::array_t<double>(
            {s1, s2},                              // Shape of the array
            {sizeof(double) * s2, sizeof(double)}, // Strides of the array (assuming it's a 1D array)
            double_ptr,                            // Pointer to the data
            py::capsule(double_ptr, [](void *double_ptr)
                        {
                            // Deleter function, called when the Python object is destroyed
                            // Do nothing here because we don't own the memory
                        }));
    }
}

void set_double_ptr(double *double_ptr, const py::array_t<double> &array, const unsigned int s1, const unsigned int s2 = 0)
{
    // Ensure the array has the correct shape and format
    auto info = array.request();
    if (info.format != py::format_descriptor<double>::format() || info.itemsize != sizeof(double))
    {
        throw std::runtime_error("Incompatible buffer format.");
    }
    if (info.ndim == 1 && info.shape[0] == 6)
    {
        // Copy data from the NumPy array to the internal buffer
        std::memcpy(double_ptr, info.ptr, sizeof(double) * s1);
    }
    else if (info.ndim == 2 && info.shape[0] == 6 && info.shape[1] == 6)
    {
        // Copy data from the NumPy array to the internal buffer
        std::memcpy(double_ptr, info.ptr, sizeof(double) * s1 * s2);
    }
    else
    {
        throw std::runtime_error("Incompatible input size.");
    }
}

class CustomArray
{
public:
    py::array_t<double> arr;
    bool &has_;
    const std::string idt;

    CustomArray(const py::array_t<double> &array, bool &has, std::string Idt)
        : arr(array), has_(has), idt(Idt){};

    double &__getitem__s(int index)
    {
        auto ptr = static_cast<double *>(this->arr.request().ptr);
        if (index < this->arr.size() && index >= 0)
        {
            return ptr[index];
        }
        else
        {
            throw std::invalid_argument("Index out of range");
        }
    }

    double &__getitem__m(const std::vector<size_t> &indices)
    {
        if (indices.size() != this->arr.ndim())
        {
            throw std::invalid_argument("Incorrect number of indices");
        }

        size_t flat_index = 0;
        size_t stride = 1;
        for (size_t i = 0; i < indices.size(); ++i)
        {
            if (indices[i] >= this->arr.shape(i))
            {
                throw std::out_of_range("Index out of range");
            }
            flat_index += indices[i] * stride;
            stride *= this->arr.shape(i);
        }

        auto ptr = static_cast<double *>(this->arr.request().ptr);

        return ptr[flat_index];
    }

    void __setitem__s(int index, double value)
    {
        auto ptr = static_cast<double *>(this->arr.request().ptr);
        if (index < this->arr.size() && index >= 0)
        {
            ptr[index] = value;
            if (this->has_ == true)
            {
                if (value == 0.0 && std::memcmp(ptr, id6, sizeof(double) * 6) == 0)
                {
                    this->has_ = false;
                }
            }
            else
            {
                if (value != 0.0)
                {
                    this->has_ = true;
                }
            }
        }
        else
        {
            throw std::invalid_argument("Index out of range");
        }
    }

    void __setitem__m(const std::vector<size_t> &indices, double value)
    {
        if (indices.size() != this->arr.ndim())
        {
            throw std::invalid_argument("Incorrect number of indices");
        }

        size_t flat_index = 0;
        size_t stride = 1;
        for (size_t i = 0; i < indices.size(); ++i)
        {
            if (indices[i] >= this->arr.shape(i))
            {
                throw std::out_of_range("Index out of range");
            }
            flat_index += indices[i] * stride;
            stride *= this->arr.shape(i);
        }

        auto ptr = static_cast<double *>(this->arr.request().ptr);

        ptr[flat_index] = value;
        if (this->has_ == true)
        {
            if (((flat_index % 7 == 0 && value == 1.0) || (flat_index % 7 != 0 && value == 0.0)) && std::memcmp(ptr, id66, sizeof(double) * 6 * 6) == 0)
            {
                this->has_ = false;
            }
        }
        else
        {
            if ((flat_index % 7 == 0 && value != 1.0) || (flat_index % 7 != 0 && value != 0.0))
            {
                this->has_ = true;
            }
        }
    }

    std::string repr() const
    {
        std::ostringstream oss;
        oss << this->idt << "(" << this->arr.str() << ")";
        return oss.str();
    }
};

void set_double_vec(std::vector<double> &double_vec, const py::array_t<double> &array)
{
    // Ensure the array has the correct shape and format
    auto info = array.request();
    if (info.ndim != 1 || info.format != py::format_descriptor<double>::format() || info.itemsize != sizeof(double))
    {
        throw std::runtime_error("Incompatible buffer format.");
    }
    if (info.shape[0] != double_vec.size())
    {
        double_vec.resize(info.shape[0]);
    }
    std::memcpy(double_vec.data(), info.ptr, sizeof(double) * info.shape[0]);
}

template <typename T>
std::vector<Pos<T>> numpyArrayToPos(const py::array_t<T> &arr)
{
    auto buffer = arr.request();
    if (buffer.ndim == 1)
    {
        if (buffer.size != 6)
        {
            throw std::runtime_error("Input array must have 6 elements");
        }

        auto ptr = static_cast<T *>(buffer.ptr);
        return {Pos<T>(ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5])};
    }
    else
    {
        if (buffer.ndim != 2 || buffer.shape[0] != 6)
        {
            throw std::runtime_error("Input array must be 2D with 6 rows");
        }

        size_t cols = buffer.shape[1];
        std::vector<Pos<T>> vecpos;
        vecpos.reserve(cols);

        auto ptr = static_cast<T *>(buffer.ptr);
        for (size_t i = 0; i < cols; ++i)
        {
            size_t offset = i * buffer.shape[0];
            vecpos.emplace_back(ptr[offset], ptr[offset + 1], ptr[offset + 2], ptr[offset + 3], ptr[offset + 4], ptr[offset + 5]);
        }

        return vecpos;
    }
}

template <typename T>
py::array_t<T> posToNumpyArray(const Pos<T> &pos)
{
    std::vector<T> data = {pos.rx, pos.px, pos.ry, pos.py, pos.de, pos.dl};
    return py::array_t<T>(data.size(), data.data());
}

template <typename T>
py::array_t<T> posToNumpyArray(const std::vector<Pos<T>> &vecpos)
{
    size_t rows = 6;             // Number of rows (attributes of Pos)
    size_t cols = vecpos.size(); // Number of columns (number of Pos objects)

    /// Create a buffer to hold the data
    std::vector<T> buffer(rows * cols);

    // Copy the data from vecpos to the buffer
    for (size_t i = 0; i < cols; ++i)
    {
        const Pos<T> &pos = vecpos[i];
        buffer[i] = pos.rx;
        buffer[i + cols] = pos.px;
        buffer[i + 2 * cols] = pos.ry;
        buffer[i + 3 * cols] = pos.py;
        buffer[i + 4 * cols] = pos.de;
        buffer[i + 5 * cols] = pos.dl;
    }

    py::array_t<T> result;
    if (cols == 1)
    {
        result = py::array_t<T>({rows}, buffer.data());
    }
    else
    {
        result = py::array_t<T>({rows, cols}, buffer.data());
    }

    return result;
}

py::array_t<double> _element_pass(const Element &e, const py::array_t<double> &array, const Accelerator &a)
{
    std::vector<Pos<double>> p = numpyArrayToPos(array);
    if (p.size() == 1)
    {
        Pos<double> p_ = p[0];
        track_elementpass(a, e, p_);
        return posToNumpyArray(p_);
    }
    else
    {
        track_elementpass(a, e, p);
        return posToNumpyArray(p);
    }
}

const std::vector<std::string> plane_dict = {"no_plane", "plane_x", "plane_y", "plane_xy"};

std::tuple<py::array_t<double>, int> _line_pass(const Accelerator &a, const py::array_t<double> &array, std::vector<unsigned int> &indices, const int &element_offset = 0)
{
    std::vector<Pos<double>> p = numpyArrayToPos(array);
    unsigned int eoffset = 0 + element_offset;
    if (p.size() == 1)
    {
        Pos<double> p_ = p[0];
        std::vector<Pos<double>> p_out;
        Plane::type lost_plane = Plane::no_plane;
        Status::type status = track_linepass(a, p_, indices, eoffset, p_out, lost_plane);
        return std::make_tuple(posToNumpyArray(p_out), 0);
    }
    else
    {
        std::vector<Pos<double>> p_out;
        std::vector<unsigned int> lost_plane;
        std::vector<bool> lost_flag;
        std::vector<int> lost_element;
        Status::type status = track_linepass(a, p, indices, eoffset, p_out, lost_plane, lost_flag, lost_element);
        return std::make_tuple(posToNumpyArray(p_out), 0);
    }
}

// template <typename T>
// py::array_t<double> _ring_pass(const Element &e, const py::array_t<double>& array, const Accelerator &a){
//     std::vector<Pos<double>> p = numpyArrayToPos(array);
//     if (p.size() == 1){
//         Pos<T> p_ = p[0];
//         track_elementpass(e, p_, a);
//         return posToNumpyArray(p_);
//     } else {
//         track_elementpass(e, p, a);
//         return posToNumpyArray(p);
//     }
// }

PYBIND11_MAKE_OPAQUE(std::vector<Element>);
PYBIND11_MODULE(trackcpp_pybind, m)
{
    m.doc() = "Trackcpp python package created with PyBind11";
    py::class_<CustomArray>(m, "CustomArray")
        .def(py::init<const py::array_t<double> &, bool &, std::string>())
        .def("__getitem__", &CustomArray::__getitem__s)
        .def("__setitem__", &CustomArray::__setitem__s)
        .def("__getitem__", &CustomArray::__getitem__m)
        .def("__setitem__", &CustomArray::__setitem__m)
        .def("__repr__", &CustomArray::repr)
    ;
    py::class_<Element>(m, "Element")
        .def(py::init<const std::string &, const double &>(),
             py::arg("fam_name") = "",
             py::arg("length") = 0.0)
        .def_readwrite("fam_name", &Element::fam_name)
        .def_readwrite("length", &Element::length)
        .def_property("pass_method", &Element::get_pass_method, &Element::set_pass_method)
        .def_readwrite("nr_steps", &Element::nr_steps)
        .def_readwrite("vchamber", &Element::vchamber)
        .def_readwrite("hmin", &Element::hmin)
        .def_readwrite("hmax", &Element::hmax)
        .def_readwrite("vmin", &Element::vmin)
        .def_readwrite("vmax", &Element::vmax)
        .def_readwrite("hkick", &Element::hkick)
        .def_readwrite("vkick", &Element::vkick)
        .def_readwrite("angle", &Element::angle)
        .def_readwrite("angle_in", &Element::angle_in)
        .def_readwrite("angle_out", &Element::angle_out)
        .def_readwrite("gap", &Element::gap)
        .def_readwrite("fint_in", &Element::fint_in)
        .def_readwrite("fint_out", &Element::fint_out)
        .def_readwrite("thin_KL", &Element::thin_KL)
        .def_readwrite("thin_SL", &Element::thin_SL)
        .def_readwrite("frequency", &Element::frequency)
        .def_readwrite("voltage", &Element::voltage)
        .def_readwrite("phase_lag", &Element::phase_lag)
        .def_property(
            "t_in",
            [](Element &element)
            {
                return CustomArray(get_double_ptr(element.t_in, 6), element.has_t_in, "T");
            },
            [](Element &element, const py::array_t<double> &array)
            {
                set_double_ptr(element.t_in, array, 6);
                if (element.has_t_in == true && std::memcmp(element.t_in, id6, sizeof(double) * 6) == 0)
                {
                    element.has_t_in = false;
                }
                else
                {
                    if (std::memcmp(element.t_in, id6, sizeof(double) * 6) != 0)
                    {
                        element.has_t_in = true;
                    }
                }
            },
            py::return_value_policy::reference_internal)
        .def_readonly("has_t_in", &Element::has_t_in)
        .def_property(
            "t_out",
            [](Element &element)
            {
                return CustomArray(get_double_ptr(element.t_out, 6), element.has_t_out, "T");
            },
            [](Element &element, const py::array_t<double> &array)
            {
                set_double_ptr(element.t_out, array, 6);
                if (element.has_t_out == true && std::memcmp(element.t_out, id6, sizeof(double) * 6) == 0)
                {
                    element.has_t_out = false;
                }
                else
                {
                    if (std::memcmp(element.t_out, id6, sizeof(double) * 6) != 0)
                    {
                        element.has_t_out = true;
                    }
                }
            },
            py::return_value_policy::reference_internal)
        .def_readonly("has_t_out", &Element::has_t_out)
        .def_property(
            "r_in",
            [](Element &element)
            {
                return CustomArray(get_double_ptr(element.r_in, 6, 6), element.has_r_in, "R");
            },
            [](Element &element, const py::array_t<double> &array)
            {
                set_double_ptr(element.r_in, array, 6, 6);
                if (element.has_r_in == true && std::memcmp(element.r_in, id66, sizeof(double) * 6 * 6) == 0)
                {
                    element.has_r_in = false;
                }
                else
                {
                    if (std::memcmp(element.r_in, id66, sizeof(double) * 6 * 6) != 0)
                    {
                        element.has_r_in = true;
                    }
                }
            },
            py::return_value_policy::reference_internal)
        .def_readonly("has_r_in", &Element::has_r_in)
        .def_property(
            "r_out",
            [](Element &element)
            {
                return CustomArray(get_double_ptr(element.r_out, 6, 6), element.has_r_out, "R");
            },
            [](Element &element, const py::array_t<double> &array)
            {
                set_double_ptr(element.r_out, array, 6, 6);
                if (element.has_r_out == true && std::memcmp(element.r_out, id66, sizeof(double) * 6 * 6) == 0)
                {
                    element.has_r_out = false;
                }
                else
                {
                    if (std::memcmp(element.r_out, id66, sizeof(double) * 6 * 6) != 0)
                    {
                        element.has_r_out = true;
                    }
                }
            },
            py::return_value_policy::reference_internal)
        .def_readonly("has_r_out", &Element::has_r_out)
        .def_property(
            "polynom_a",
            [](const Element &e)
            {
                return get_double_ptr(e.polynom_a.data(), e.polynom_a.size());
            },
            [](Element &element, const py::array_t<double> &array)
            {
                set_double_vec(element.polynom_a, array);
            })
        .def_property(
            "polynom_b",
            [](const Element &e)
            {
                return get_double_ptr(e.polynom_b.data(), e.polynom_b.size());
            },
            [](Element &element, const py::array_t<double> &array)
            {
                set_double_vec(element.polynom_b, array);
            })
        .def_property(
            "K",
            [](const Element &e){
                return e.polynom_b[1];
            },
            [](Element &e, const double value){
                e.polynom_b[1] = value;
            })
        .def_property(
            "KL",
            [](const Element &e){
                return e.polynom_b[1]*e.length;
            },
            [](Element &e, const double value){
                e.polynom_b[1] = value/e.length;
            })
        .def_property(
            "KxL",
            [](const Element &e){
                return -e.matrix66[1][0];
            },
            [](Element &e, const double value){
                e.matrix66[1][0] = -value;
            })
        .def_property(
            "KyL",
            [](const Element &e){
                return -e.matrix66[3][2];
            },
            [](Element &e, const double value){
                e.matrix66[3][2] = -value;
            })
        .def_property(
            "Ks",
            [](const Element &e){
                return -e.polynom_a[1];
            },
            [](Element &e, const double value){
                e.polynom_a[1] = -value;
            })
        .def_property(
            "KsL",
            [](const Element &e){
                return -e.polynom_a[1]*e.length;
            },
            [](Element &e, const double value){
                e.polynom_a[1] = -value/e.length;
            })
        .def_property(
            "KsxL",
            [](const Element &e){
                return -e.matrix66[1][2];
            },
            [](Element &e, const double value){
                e.matrix66[1][2] = -value;
            })
        .def_property(
            "KsyL",
            [](const Element &e){
                return -e.matrix66[3][0];
            },
            [](Element &e, const double value){
                e.matrix66[3][0] = -value;
            })
        .def_property(
            "S",
            [](const Element &e){
                return e.polynom_b[2];
            },
            [](Element &e, const double value){
                e.polynom_b[2] = value;
            })
        .def_property(
            "SL",
            [](const Element &e){
                return e.polynom_b[2]*e.length;
            },
            [](Element &e, const double value){
                e.polynom_b[2] = value/e.length;
            })
        .def_property(
            "hkick_polynom",
            [](const Element &e){
                return e.polynom_b[0] * (-e.length);
            },
            [](Element &e, const double value){
                e.polynom_b[0] = -value/e.length;
            }
        )
        .def_property(
            "vkick_polynom",
            [](const Element &e){
                return e.polynom_b[0] * e.length;
            },
            [](Element &e, const double value){
                e.polynom_b[0] = value/e.length;
            }
        )
        .def("__eq__",
            [](const Element &self, const Element &other){
                return (self == other);
            }
        )
        .def("__eq__",
            [](const Element &self, py::object other){
                return false;
            }
        )
        .def("__repr__",
            [](const Element &e){
                std::ostringstream rst;
                rst << "fam_name: " << e.fam_name;
                return rst.str();
            }
        )
        .def("__str__",
            [](const Element &e){
                std::ostringstream rst;
                rst << e;
                return rst.str();
            }
        )
        // .def_property("matrix66",
        // [](Element& e){
        //     return CustomMatrix(e.matrix66);
        // },
        // [](Element& e, const py::array_t<double>& mat){}
        // )
        ;

    // Elements
    m.def("marker", &Element::marker);
    m.def("bpm", &Element::bpm);
    m.def("hcorrector", &Element::hcorrector);
    m.def("vcorrector", &Element::vcorrector);
    m.def("corrector", &Element::corrector);
    m.def("drift", &Element::drift);
    m.def("drift_g2l", &Element::drift_g2l);
    m.def("matrix", &Element::matrix);
    m.def("rbend", &Element::rbend);
    m.def("quadrupole", &Element::quadrupole);
    m.def("sextupole", &Element::sextupole);
    m.def("rfcavity", &Element::rfcavity);
    m.def("kickmap", &Element::kickmap);

    py::bind_vector<std::vector<Element>>(m, "CppElementVector");

    // Accelerator
    py::class_<Accelerator>(m, "Accelerator")
        .def(py::init<const double &>(),
             py::arg("energy") = -1)
        .def_property(
            "energy",
            [](const Accelerator &a)
            {
                return a.energy;
            },
            [](Accelerator &a, const double &energy)
            {
                if (energy > electron_rest_energy_eV)
                {
                    a.energy = energy;
                }
            })
        .def_property(
            "cavity_on",
            [](const Accelerator &a)
            {
                return bool(a.cavity_on);
            },
            [](Accelerator &a, const bool value)
            {
                a.cavity_on = value;
            })
        .def_property(
            "radiation_on",
            [](const Accelerator &a)
            {
                return a.radiation_on;
            },
            [](Accelerator &a, const int value)
            {
                if (value <= 2 && 0 <= value)
                {
                    a.radiation_on = value;
                }
            })
        .def_property(
            "vchamber_on",
            [](const Accelerator &a)
            {
                return a.vchamber_on;
            },
            [](Accelerator &a, const bool value)
            {
                a.vchamber_on = value;
            })
        .def_readwrite("harmonic_number", &Accelerator::harmonic_number)
        .def_readwrite("lattice_version", &Accelerator::lattice_version)
        .def_property_readonly("length", [](const Accelerator &a)
                               { return a.get_length(); })
        .def("append",
             [](Accelerator &a, const Element &e)
             {
                 a.lattice.push_back(Element(e));
             })
        .def("insert",
             [](Accelerator &a, const int index, const Element &e)
             {
                 a.lattice.insert(a.lattice.begin() + index, Element(e));
             })
        .def(
            "pop",
            [](Accelerator &a)
            {
                Element to_pop = Element(a.lattice.back());
                a.lattice.pop_back();
                return to_pop;
            },
            py::return_value_policy::reference_internal)
        .def(
            "__getitem__", [](const Accelerator &a, const int index)
            {
            size_t s = a.lattice.size();
            if (s > 0 && index < s) {
                return a.lattice[index];
            } else {
                throw std::out_of_range("Index out of Lattice bounds.");
            } },
            py::return_value_policy::reference_internal)
        .def("__len__", [](const Accelerator &a)
             { return a.lattice.size(); });

    // tracking
    m.def("element_pass", &_element_pass, ".");
    m.def("line_pass", &_line_pass, ".");
    // m.def("ring_pass", &_ring_pass, ".");
}
