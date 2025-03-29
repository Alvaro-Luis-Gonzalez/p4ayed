// AUTOR: 
// FECHA: 
// EMAIL: 
// VERSION: 2.0
// ASIGNATURA: Algoritmos y Estructuras de Datos
// PRÁCTICA Nº: 4
// ESTILO: Google C++ Style Guide
// COMENTARIOS:
// 

#ifndef SLLPOLYNOMIAL_H_
#define SLLPOLYNOMIAL_H_

#include <iostream>
#include <math.h>  // fabs, pow

#include "pair_t.h"
#include "sll_t.h"
#include "vector_t.h"

#define EPS 1.0e-6

typedef pair_t<double> pair_double_t;  // Campo data_ de SllPolynomial
typedef sll_node_t<pair_double_t> SllPolyNode;  // Nodos de SllPolynomial

// Clase para polinomios basados en listas simples de pares
class SllPolynomial : public sll_t<pair_double_t> {
 public:
  // constructores
  SllPolynomial(void) : sll_t() {};
  SllPolynomial(const vector_t<double>&, const double = EPS);

  // destructor
  ~SllPolynomial() {};

  // E/S
  void Write(std::ostream& = std::cout) const;
  
  // operaciones
  double Eval(const double) const;
  bool IsEqual(const SllPolynomial&, const double = EPS) const;
  void Sum(const SllPolynomial&, SllPolynomial&, const double = EPS);
};


bool IsNotZero(const double val, const double eps = EPS) {
  return fabs(val) > eps;
}

// FASE II
// constructor
SllPolynomial::SllPolynomial(const vector_t<double>& v, const double eps) {
  for (int i = 0; i < v.get_size(); ++i) {
    double coeficiente = v[i];
    // Verificamos si el coeficiente es distinto de cero (considerando la precisión)
    if (IsNotZero(coeficiente, eps)) {
      // Creamos un nuevo nodo para el término del polinomio
      SllPolyNode* nuevo_nodo = new SllPolyNode(pair_double_t(coeficiente, i));
      // Insertamos el nuevo nodo al final de la lista
      // Para insertar al final, hay que recorrer la lista. 
      // Si la lista está vacía, el nuevo nodo será la cabeza
      if (this->empty()) {
        this->push_front(nuevo_nodo); // Usamos push_front para insertar al inicio (es más eficiente)
      } else {
        SllPolyNode* aux = this->get_head();
        while (aux->get_next() != NULL) {
          aux = aux->get_next();
        }
        this->insert_after(aux, nuevo_nodo);
      }
    }
  }
}
  
// E/S
void SllPolynomial::Write(std::ostream& os) const {
  os << "[ ";
  bool first{true};
  SllPolyNode* aux{get_head()};
  while (aux != NULL) {
    int inx{aux->get_data().get_inx()};
    double val{aux->get_data().get_val()};
    if (val > 0)
      os << (!first ? " + " : "") << val;
    else
      os << (!first ? " - " : "-") << fabs(val);
    os << (inx > 1 ? " x^" : (inx == 1) ? " x" : "");
    if (inx > 1)
      os << inx;
    first = false;
    aux = aux->get_next();
  }
  os << " ]" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SllPolynomial& p) {
  p.Write(os);
  return os;
}


// Operaciones con polinomios

// FASE III
// Evaluación de un polinomio representado por lista simple
double SllPolynomial::Eval(const double x) const {
  double result{0.0};
  SllPolyNode* aux = this->get_head();

  while (aux != NULL) {
    double coeficiente = aux->get_data().get_val();
    int exponente = aux->get_data().get_inx();
    result += coeficiente * pow(x, exponente);
    aux = aux->get_next();
  }

  return result;
}

// Comparación si son iguales dos polinomios representados por listas simples
bool SllPolynomial::IsEqual(const SllPolynomial& sllpol, const double eps) const {
  bool differents = false;

  SllPolyNode* aux1 = this->get_head();    // Primer polinomio
  SllPolyNode* aux2 = sllpol.get_head();   // Segundo polinomio

  while (aux1 != NULL && aux2 != NULL && !differents) {
    // Comparar exponentes
    if (aux1->get_data().get_inx() != aux2->get_data().get_inx()) {
      differents = true;
    }

    // Comparar coeficientes, considerando valores pequeños como cero
    double coef1 = aux1->get_data().get_val();
    double coef2 = aux2->get_data().get_val();

    if (!IsNotZero(coef1, eps)) coef1 = 0.0;
    if (!IsNotZero(coef2, eps)) coef2 = 0.0;

    if (fabs(coef1 - coef2) > eps) {
      differents = true;
    }

    // Avanzamos en ambas listas
    aux1 = aux1->get_next();
    aux2 = aux2->get_next();
  }

  // Si una lista tiene más términos que la otra, son diferentes
  if ((aux1 != NULL || aux2 != NULL)) {
    differents = true;
  }

  return !differents;
}


// FASE IV
// Generar nuevo polinomio suma del polinomio invocante mas otro polinomio
void SllPolynomial::Sum(const SllPolynomial& sllpol,
			SllPolynomial& sllpolsum,
			const double eps) {

SllPolyNode* aux_this = this->get_head();
SllPolyNode* aux_other = sllpol.get_head();

// Recorremos ambas listas simultáneamente
while (aux_this != NULL || aux_other != NULL) {
  double coeficiente_suma = 0.0;
  int exponente;

  // Caso 1: Solo quedan términos en el primer polinomio
  if (aux_this != NULL && aux_other == NULL) {
    exponente = aux_this->get_data().get_inx();
    coeficiente_suma = aux_this->get_data().get_val();
    aux_this = aux_this->get_next();
  }
  // Caso 2: Solo quedan términos en el segundo polinomio
  else if (aux_this == NULL && aux_other != NULL) {
    exponente = aux_other->get_data().get_inx();
    coeficiente_suma = aux_other->get_data().get_val();
    aux_other = aux_other->get_next();
  }
  // Caso 3: Términos con el mismo exponente en ambos polinomios
  else if (aux_this != NULL && aux_other != NULL &&
           aux_this->get_data().get_inx() == aux_other->get_data().get_inx()) {
    exponente = aux_this->get_data().get_inx();
    coeficiente_suma = aux_this->get_data().get_val() +
                       aux_other->get_data().get_val();
    aux_this = aux_this->get_next();
    aux_other = aux_other->get_next();
  }
  // Caso 4: Término con menor exponente en el primer polinomio
  else if (aux_this != NULL && aux_other != NULL &&
           aux_this->get_data().get_inx() < aux_other->get_data().get_inx()) {
    exponente = aux_this->get_data().get_inx();
    coeficiente_suma = aux_this->get_data().get_val();
    aux_this = aux_this->get_next();
  }
  // Caso 5: Término con menor exponente en el segundo polinomio
  else {
    exponente = aux_other->get_data().get_inx();
    coeficiente_suma = aux_other->get_data().get_val();
    aux_other = aux_other->get_next();
  }

  // Insertamos el término en el polinomio suma (si el coeficiente no es cero)
  if (IsNotZero(coeficiente_suma, eps)) {
    SllPolyNode* nuevo_nodo =
        new SllPolyNode(pair_double_t(coeficiente_suma, exponente));
    // Insertar al final (mismo criterio que en el constructor para simplificar)
    if (sllpolsum.empty()) {
      sllpolsum.push_front(nuevo_nodo);
    } else {
      SllPolyNode* aux = sllpolsum.get_head();
      while (aux->get_next() != NULL) {
        aux = aux->get_next();
      }
      sllpolsum.insert_after(aux, nuevo_nodo);
    }
  }
}
}       

#endif  // SLLPOLYNOMIAL_H_
