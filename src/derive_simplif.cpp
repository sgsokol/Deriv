#include <Rcpp.h>
#include <string>
#include <map>
#include <vector>
#include <cmath>

using namespace Rcpp;
// Compatibility macros for R < 4.5.0
#if R_VERSION < R_Version(4, 5, 0)
#  define R_ClosureFormals(x) FORMALS(x)
#  define R_ClosureBody(x)    BODY(x)
#  define R_ClosureEnv(x)     CLOENV(x)
#endif

void debug_print(const char* label, SEXP expr) {
    Rcpp::Rcout << "[DEBUG] " << label << ": ";
    Rf_PrintValue(expr);
}

// Returns true if expr is a unary minus call, e.g. (-x)
inline bool is_unary_minus(SEXP expr) {
    return TYPEOF(expr) == LANGSXP && 
           Rf_length(expr) == 2 && 
           CAR(expr) == Rf_install("-");
}

// Helper to substitute formal parameter names with actual calling arguments inside a drule template
SEXP substitute_expr(SEXP expr, Rcpp::List subst_map) {
    if (expr == R_NilValue) {
        return R_NilValue;
    }

    // 1. Symbol Substitution
    if (TYPEOF(expr) == SYMSXP) {
        const char* name_str = CHAR(PRINTNAME(expr));
        if (subst_map.containsElementNamed(name_str)) {
            return subst_map[name_str];
        }
        return expr;
    }

    // 2. Language / Call Tree Substitution
    if (TYPEOF(expr) == LANGSXP) {
        SEXP res = PROTECT(Rf_shallow_duplicate(expr));
        SEXP curr_old = expr;
        SEXP curr_new = res;
        while (curr_old != R_NilValue) {
            SETCAR(curr_new, substitute_expr(CAR(curr_old), subst_map));
            SET_TAG(curr_new, TAG(curr_old));
            curr_old = CDR(curr_old);
            curr_new = CDR(curr_new);
        }
        UNPROTECT(1);
        return res;
    }
    return expr;
}

// Helper functions for AST inspection and formatting
bool is_numeric_val(SEXP x, double val) {
    if (TYPEOF(x) == REALSXP && Rf_length(x) == 1 && REAL(x)[0] == val) return true;
    if (TYPEOF(x) == INTSXP && Rf_length(x) == 1 && INTEGER(x)[0] == val) return true;
    return false;
}

bool get_numeric_val(SEXP x, double& val) {
    if (TYPEOF(x) == REALSXP && Rf_length(x) == 1) {
        val = REAL(x)[0];
        return true;
    }
    if (TYPEOF(x) == INTSXP && Rf_length(x) == 1) {
        val = (double)INTEGER(x)[0];
        return true;
    }
    return false;
}

std::string deparse_node(SEXP expr) {
    if (TYPEOF(expr) == SYMSXP) {
        return CHAR(PRINTNAME(expr));
    }
    if (TYPEOF(expr) == LANGSXP) {
        SEXP fun = CAR(expr);
        if (TYPEOF(fun) == SYMSXP) {
            std::string op = CHAR(PRINTNAME(fun));
            int len = Rf_length(expr);
            if ((op == "[" || op == "[[" || op == "$") && len >= 3) {
                std::string obj = deparse_node(CADR(expr));
                std::string idx = deparse_node(CADDR(expr));
                return obj + "[" + idx + "]";
            }
        }
    }
    if (TYPEOF(expr) == REALSXP && Rf_length(expr) == 1) {
        double v = REAL(expr)[0];
        if (v == std::floor(v)) return std::to_string((int)v);
        return std::to_string(v);
    }
    if (TYPEOF(expr) == INTSXP && Rf_length(expr) == 1) {
        return std::to_string(INTEGER(expr)[0]);
    }
    if (TYPEOF(expr) == STRSXP && Rf_length(expr) == 1) {
        return std::string(CHAR(STRING_ELT(expr, 0)));
    }
    return "";
}

// ---------------------------------------------------------
// Simplify Pass
// ---------------------------------------------------------

//' Simplify Symbolic Expressions
//'
//' @param expr An R language object, symbol, or numeric constant to simplify.
//' @return A simplified R expression or scalar numeric value.
//' @export
// [[Rcpp::export]]
SEXP Simplify_cpp(SEXP expr) {
    if (TYPEOF(expr) != LANGSXP) return expr;

    SEXP fun = CAR(expr);
    if (TYPEOF(fun) != SYMSXP) return expr; 
    
    std::string op = CHAR(PRINTNAME(fun));
    int len = Rf_length(expr);
    
    if (len == 3) {
        // Evaluate and protect left and right children
        SEXP left = PROTECT(Simplify_cpp(CADR(expr)));
        SEXP right = PROTECT(Simplify_cpp(CADDR(expr)));
        
        double val_l, val_r;
        if (get_numeric_val(left, val_l) && get_numeric_val(right, val_r)) {
            SEXP res = R_NilValue;
            if (op == "+") res = Rcpp::wrap(val_l + val_r);
            else if (op == "-") res = Rcpp::wrap(val_l - val_r);
            else if (op == "*") res = Rcpp::wrap(val_l * val_r);
            else if (op == "/" && val_r != 0.0) res = Rcpp::wrap(val_l / val_r);
            else if (op == "^") res = Rcpp::wrap(std::pow(val_l, val_r));
            
            if (res != R_NilValue) {
                UNPROTECT(2); // left, right
                return res;
            }
        }
        
        if (op == "+") {
            if (is_numeric_val(left, 0.0))  { UNPROTECT(2); return right; }
            if (is_numeric_val(right, 0.0)) { UNPROTECT(2); return left; }

            // x + (-y)  ==>  x - y
            if (is_unary_minus(right)) {
                SEXP inner_y = CADR(right);
                SEXP new_sub = PROTECT(Rcpp::Language("-", left, inner_y));
                SEXP res = Simplify_cpp(new_sub);
                UNPROTECT(3); // left, right, new_sub
                return res;
            }

            // (-y) + x  ==>  x - y
            if (is_unary_minus(left)) {
                SEXP inner_y = CADR(left);
                SEXP new_sub = PROTECT(Rcpp::Language("-", right, inner_y));
                SEXP res = Simplify_cpp(new_sub);
                UNPROTECT(3); // left, right, new_sub
                return res;
            }

        } else if (op == "*") {
            if (is_numeric_val(left, 0.0) || is_numeric_val(right, 0.0)) {
                UNPROTECT(2);
                return Rcpp::wrap(0.0);
            }
            if (is_numeric_val(left, 1.0))  { UNPROTECT(2); return right; }
            if (is_numeric_val(right, 1.0)) { UNPROTECT(2); return left; }
            
            if (is_numeric_val(right, -1.0)) {
                SEXP neg_call = PROTECT(Rcpp::Language("-", left));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(3); // left, right, neg_call
                return res;
            }
            if (is_numeric_val(left, -1.0)) {
                SEXP neg_call = PROTECT(Rcpp::Language("-", right));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(3); // left, right, neg_call
                return res;
            }

            // (-a) * b ==> -(a * b)
            if (is_unary_minus(left)) {
                SEXP inner_a = CADR(left);
                SEXP prod = PROTECT(Rcpp::Language("*", inner_a, right));
                SEXP neg_call = PROTECT(Rcpp::Language("-", Simplify_cpp(prod)));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(4); // left, right, prod, neg_call
                return res;
            }

            // a * (-b) ==> -(a * b)
            if (is_unary_minus(right)) {
                SEXP inner_b = CADR(right);
                SEXP prod = PROTECT(Rcpp::Language("*", left, inner_b));
                SEXP neg_call = PROTECT(Rcpp::Language("-", Simplify_cpp(prod)));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(4); // left, right, prod, neg_call
                return res;
            }

        } else if (op == "-") {
            if (is_numeric_val(right, 0.0)) { UNPROTECT(2); return left; }
            if (is_numeric_val(left, 0.0)) {
                SEXP neg_call = PROTECT(Rcpp::Language("-", right));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(3); // left, right, neg_call
                return res;
            }

        } else if (op == "/") {
            if (is_numeric_val(left, 0.0))  { UNPROTECT(2); return Rcpp::wrap(0.0); }
            if (is_numeric_val(right, 1.0)) { UNPROTECT(2); return left; }

            // (-a) / b ==> -(a / b)
            if (is_unary_minus(left)) {
                SEXP inner_a = CADR(left);
                SEXP div = PROTECT(Rcpp::Language("/", inner_a, right));
                SEXP neg_call = PROTECT(Rcpp::Language("-", Simplify_cpp(div)));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(4); // left, right, div, neg_call
                return res;
            }

            // a / (-b) ==> -(a / b)
            if (is_unary_minus(right)) {
                SEXP inner_b = CADR(right);
                SEXP div = PROTECT(Rcpp::Language("/", left, inner_b));
                SEXP neg_call = PROTECT(Rcpp::Language("-", Simplify_cpp(div)));
                SEXP res = Simplify_cpp(neg_call);
                UNPROTECT(4); // left, right, div, neg_call
                return res;
            }

        } else if (op == "^") {
            if (is_numeric_val(right, 0.0)) { UNPROTECT(2); return Rcpp::wrap(1.0); }
            if (is_numeric_val(right, 1.0)) { UNPROTECT(2); return left; }
            if (is_numeric_val(left, 0.0))  { UNPROTECT(2); return Rcpp::wrap(0.0); }
            if (is_numeric_val(left, 1.0))  { UNPROTECT(2); return Rcpp::wrap(1.0); }
        }

        // Reconstruct call node if unchanged
        SEXP res = Rcpp::Language(op.c_str(), left, right);
        UNPROTECT(2); // left, right
        return res;
        
    } else if (len == 2) {
        SEXP arg = PROTECT(Simplify_cpp(CADR(expr)));
        
        if (op == "+" || op == "(") {
            UNPROTECT(1);
            return arg;
        }
        
        if (op == "-") {
            double val_arg;
            if (get_numeric_val(arg, val_arg)) {
                UNPROTECT(1);
                return Rcpp::wrap(-val_arg);
            }
            
            // Double negation: -(-y) ==> y
            if (is_unary_minus(arg)) {
                SEXP inner_y = CADR(arg);
                SEXP res = Simplify_cpp(inner_y);
                UNPROTECT(1); // arg
                return res;
            }
        }

        SEXP res = Rcpp::Language(op.c_str(), arg);
        UNPROTECT(1); // arg
        return res;
    }
    
    return expr;
}

SEXP deriv_cpp_internal(SEXP expr, const std::string& target, Rcpp::Environment drule, Rcpp::Environment env) {
    // 1. LEAF NODES
    if (TYPEOF(expr) == REALSXP || TYPEOF(expr) == INTSXP) {
        return Rf_ScalarReal(0.0);
    }
    
    if (TYPEOF(expr) == SYMSXP) {
        std::string var_name = CHAR(PRINTNAME(expr));
        return Rf_ScalarReal(var_name == target ? 1.0 : 0.0);
    }

    // 2. CALL NODES (LANGSXP)
    if (TYPEOF(expr) == LANGSXP) {
        SEXP fun = CAR(expr);
        if (TYPEOF(fun) != SYMSXP) {
            Rcpp::stop("Function call head is not a valid symbol");
        }
        
        std::string op = CHAR(PRINTNAME(fun));
        int len = Rf_length(expr);

        // --- ACCESSORS & INDEXING (`[`, `[[`, `$`) ---
        if (op == "[" || op == "[[" || op == "$") {
            std::string expr_str = deparse_node(expr);
            if (!expr_str.empty() && expr_str == target) {
                return Rf_ScalarReal(1.0);
            }
            return Rf_ScalarReal(0.0);
        }

        // --- A. SPECIAL OPERATOR CASES ---
        if (op == "(" && len == 2) {
            return deriv_cpp_internal(CADR(expr), target, drule, env);
        }

        if (len == 2) {
            if (op == "+") {
                return deriv_cpp_internal(CADR(expr), target, drule, env);
            }
            if (op == "-") {
                SEXP d_arg = PROTECT(deriv_cpp_internal(CADR(expr), target, drule, env));
                SEXP res = Simplify_cpp(Rcpp::Language("-", d_arg));
                UNPROTECT(1);
                return res;
            }
        }

        if (len == 3) {
            if (op == "+" || op == "-") {
                SEXP d_left  = PROTECT(deriv_cpp_internal(CADR(expr), target, drule, env));
                SEXP d_right = PROTECT(deriv_cpp_internal(CADDR(expr), target, drule, env));
                
                SEXP res = Simplify_cpp(Rcpp::Language(op.c_str(), d_left, d_right));
                UNPROTECT(2);
                return res;
            }
        }

        // --- B. LOOKUP IN DRULE TABLE ---
        //debug_print("fun=", fun);
        SEXP rules = drule.get(CHAR(PRINTNAME(fun)));
        if (rules != R_NilValue && TYPEOF(rules) == VECSXP) {
            
            SEXP arg_names = Rf_getAttrib(rules, R_NamesSymbol);
            int num_formals = Rf_length(rules);
            
            SEXP subst_env = PROTECT(R_NewEnv(R_NilValue, 0, 0));
            SEXP curr_args = CDR(expr);
            int actual_len = Rf_length(curr_args);
            
            for (int i = 0; i < num_formals; i++) {
                SEXP formal_sym = Rf_install(CHAR(STRING_ELT(arg_names, i)));
                SEXP actual_val = R_NilValue;
                
                if (curr_args != R_NilValue) {
                    actual_val = CAR(curr_args);
                    curr_args = CDR(curr_args);
                } else {
                    actual_val = formal_sym; 
                }
                
                Rf_defineVar(formal_sym, actual_val, subst_env);
            }
            
            SEXP total_deriv = R_NilValue;
            curr_args = CDR(expr);
            
            for (int i = 0; i < actual_len && i < num_formals; i++) {
                SEXP rule_template = VECTOR_ELT(rules, i);
                
                if (rule_template == R_NilValue || curr_args == R_NilValue) {
                    if (curr_args != R_NilValue) curr_args = CDR(curr_args);
                    continue;
                }
                
                SEXP actual_arg = CAR(curr_args);
                SEXP d_arg = PROTECT(deriv_cpp_internal(actual_arg, target, drule, env));
                
                if (is_numeric_val(d_arg, 0.0)) {
                    UNPROTECT(1);
                    curr_args = CDR(curr_args);
                    continue;
                }
                //debug_print("subst_env", Rcpp::as<Rcpp::List>(subst_env));
                SEXP df_darg = PROTECT(substitute_expr(rule_template, subst_env));
                SEXP term = PROTECT(Rcpp::Language("*", df_darg, d_arg));
                
                if (total_deriv == R_NilValue) {
                    total_deriv = term;
                    UNPROTECT(2);
                } else {
                    total_deriv = PROTECT(Rcpp::Language("+", total_deriv, term));
                    UNPROTECT(3);
                }
                
                UNPROTECT(1); // d_arg
                curr_args = CDR(curr_args);
            }
            
            UNPROTECT(1); // subst_env
            
            if (total_deriv == R_NilValue) {
                return Rf_ScalarReal(0.0);
            }
            
            return Simplify_cpp(total_deriv);
        }
        
        // --- C. LOOKUP USER-DEFINED FUNCTION BODY IN ENV ---
        //debug_print("env", Rcpp::as<Rcpp::List>(env));
        SEXP user_fun = env.find(op);

        if (user_fun != R_UnboundValue && TYPEOF(user_fun) == CLOSXP) {
            SEXP formals = R_ClosureFormals(user_fun);
            SEXP body    = R_ClosureExpr(user_fun);

            SEXP actual_args = CDR(expr); // Skip function symbol
            
            // Create named Rcpp::List for parameter mapping
            Rcpp::List subst_map;

            while (formals != R_NilValue && actual_args != R_NilValue) {
                SEXP formal_sym = TAG(formals);

                // Fallback if formal parameter symbol is in CAR(formals)
                if (formal_sym == R_NilValue || TYPEOF(formal_sym) != SYMSXP) {
                    if (TYPEOF(CAR(formals)) == SYMSXP) {
                        formal_sym = CAR(formals);
                    }
                }
                if (formal_sym != R_NilValue && TYPEOF(formal_sym) == SYMSXP) {
                    subst_map[CHAR(PRINTNAME(formal_sym))] = CAR(actual_args);
                }
                formals = CDR(formals);
                actual_args = CDR(actual_args);
            }

            // Substitute body using the List mapping
            SEXP expanded_body = PROTECT(substitute_expr(body, subst_map));
            SEXP deriv_res = deriv_cpp_internal(expanded_body, target, drule, env);

            UNPROTECT(1); // expanded_body
            return Simplify_cpp(deriv_res);
        }

        // --- D. UNKNOWN FUNCTION ERROR ---
        Rcpp::stop("Function '%s' has no derivative rule in 'drule' nor a function definition in 'env'", op.c_str());
    }

    Rcpp::stop("Unsupported expression type encountered during differentiation");
}

//' Symbolic Differentiation (Supports Vectors of Variables)
//'
//' Computes the symbolic derivative of an R expression with respect to one
//' or more target variables. If a single target is provided, it returns a 
//' scalar expression. If multiple targets are provided, it returns a named c() call.
//' 
//' This version is much faster for simple expressions but has less features than its R counterpart Deriv().
//'
//' @param expr An R language object, symbol, or numeric constant to differentiate.
//' @param x A character vector, symbol, or named vector of targets.
//' @param env An environment to resolve function calls in (defaults to parent.frame()).
//' @param drule An environment containing derivative rules (cf. Details in \code{\link{Deriv}} for syntax rules).
//' @return A simplified R expression (or a named c() call for multiple inputs).
//' @seealso \code{\link{Deriv}}
//' @export
// [[Rcpp::export(signature = {expr, x, env=parent.frame(), drule=Deriv::drule})]]
SEXP Deriv_cpp(SEXP expr, SEXP x, Rcpp::Environment env, SEXP drule) {
    int n = 0;
    bool is_strsxp = (TYPEOF(x) == STRSXP);
    
    if (is_strsxp) {
        n = Rf_length(x);
    } else if (TYPEOF(x) == SYMSXP) {
        n = 1;
    } else {
        stop("Invalid differentiation variable format. Must be string, symbol, or named vector.");
    }

    if (n == 0) {
        stop("Differentiation variable vector cannot be empty.");
    }

    SEXP names = is_strsxp ? Rf_getAttrib(x, R_NamesSymbol) : R_NilValue;

    if (n == 1) {
        std::string target_string = "";
        if (is_strsxp) {
            std::string val = CHAR(STRING_ELT(x, 0));
            bool has_name = (names != R_NilValue && Rf_length(names) > 0 && CHAR(STRING_ELT(names, 0))[0] != '\0');
            if (has_name) {
                target_string = std::string(CHAR(STRING_ELT(names, 0))) + "[" + val + "]";
            } else {
                target_string = val;
            }
        } else {
            target_string = CHAR(PRINTNAME(x));
        }
        
        SEXP raw_deriv = deriv_cpp_internal(expr, target_string, drule, env);
        return Simplify_cpp(raw_deriv);
    }

    SEXP call = PROTECT(Rf_allocLang(n + 1));
    SETCAR(call, Rf_install("c"));
    SEXP curr_cell = CDR(call);

    for (int i = 0; i < n; ++i) {
        std::string target_string = "";
        std::string output_name = "";

        if (is_strsxp) {
            std::string val = CHAR(STRING_ELT(x, i));
            bool has_name = (names != R_NilValue && Rf_length(names) > i && CHAR(STRING_ELT(names, i))[0] != '\0');
            
            if (has_name) {
                std::string name = CHAR(STRING_ELT(names, i));
                target_string = name + "[" + val + "]";
                output_name = name + "_" + val;
            } else {
                target_string = val;
                output_name = val;
            }
        } else {
            target_string = CHAR(PRINTNAME(x));
            output_name = target_string;
        }

        SEXP raw_deriv = deriv_cpp_internal(expr, target_string, drule, env);
        SEXP simplified_deriv = Simplify_cpp(raw_deriv);

        SETCAR(curr_cell, simplified_deriv);
        if (!output_name.empty()) {
            SET_TAG(curr_cell, Rf_install(output_name.c_str()));
        }
        curr_cell = CDR(curr_cell);
    }

    UNPROTECT(1);
    return call;
}
