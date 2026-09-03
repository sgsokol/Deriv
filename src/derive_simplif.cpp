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
SEXP substitute_expr(SEXP expr, SEXP env) {
    if (expr == R_NilValue || expr == R_UnboundValue) return expr;

    switch (TYPEOF(expr)) {
        case SYMSXP: {
            SEXP val = R_getVarEx(expr, env, TRUE, R_UnboundValue);
            if (val != R_UnboundValue) {
                return val; // Substitute the variable/default value
            }
            return expr; // Leave unbound symbols as-is
        }
        case LANGSXP: {
            // Reconstruct language call recursively
            SEXP head = CAR(expr);
            SEXP new_call = PROTECT(Rf_lcons(substitute_expr(head, env), R_NilValue));
            SEXP tail_out = new_call;
            
            SEXP tail_in = CDR(expr);
            while (tail_in != R_NilValue) {
                SEXP evaluated_item = substitute_expr(CAR(tail_in), env);
                SETCDR(tail_out, Rf_lcons(evaluated_item, R_NilValue));
                
                // Preserve argument tags (names) if present
                if (TAG(tail_in) != R_NilValue) {
                    SET_TAG(tail_out, TAG(tail_in));
                }
                
                tail_out = CDR(tail_out);
                tail_in = CDR(tail_in);
            }
            UNPROTECT(1);
            return new_call;
        }
        // Pass literals through safely (including LGLSXP for TRUE/FALSE defaults)
        case REALSXP:
        case INTSXP:
        case LGLSXP:
        case STRSXP:
            return expr;
            
        default:
            //Rcout << "type: " << Rf_type2char(TYPEOF(expr)) << std::endl;
            Rf_error("Unsupported expression type encountered during differentiation");
            return R_NilValue;
    }
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

inline bool is_call_to(SEXP expr, const char* name) {
    if (TYPEOF(expr) != LANGSXP) return false;
    return TYPEOF(CAR(expr)) == SYMSXP && strcmp(CHAR(PRINTNAME(CAR(expr))), name) == 0;
}

bool nodes_equal(SEXP a, SEXP b) {
    if (a == b) return true;
    if (TYPEOF(a) != TYPEOF(b)) return false;

    switch (TYPEOF(a)) {
        case REALSXP:
            return Rf_asReal(a) == Rf_asReal(b);
        case INTSXP:
            return Rf_asInteger(a) == Rf_asInteger(b);
        case SYMSXP:
            return std::string(CHAR(PRINTNAME(a))) == std::string(CHAR(PRINTNAME(b)));
        case LANGSXP:
        case LISTSXP: {
            if (Rf_length(a) != Rf_length(b)) return false;
            if (!nodes_equal(CAR(a), CAR(b))) return false;
            if (!nodes_equal(TAG(a), TAG(b))) return false;
            return nodes_equal(CDR(a), CDR(b));
        }
        default:
            return false;
    }
}

SEXP match_call_args(SEXP formals_pairlist, SEXP actual_args_pairlist) {
    SEXP subst_env = PROTECT(R_NewEnv(R_EmptyEnv, 0, 0));

    // Step 1: Bind defaults first
    SEXP f = formals_pairlist;
    while (f != R_NilValue) {
        SEXP sym = TAG(f);
        if (TYPEOF(sym) == SYMSXP) {
            SEXP default_val = CAR(f);
            if (default_val != R_MissingArg) {
                Rf_defineVar(sym, default_val, subst_env);
            }
        }
        f = CDR(f);
    }

    // Step 2: Bind supplied arguments (handles both positional and named)
    SEXP a = actual_args_pairlist;
    SEXP f_pos = formals_pairlist;

    while (a != R_NilValue) {
        SEXP arg_val = CAR(a);
        SEXP arg_tag = TAG(a);

        if (arg_tag != R_NilValue && TYPEOF(arg_tag) == SYMSXP) {
            // Named argument: e.g., log(x, base = 10)
            Rf_defineVar(arg_tag, arg_val, subst_env);
        } else if (f_pos != R_NilValue) {
            // Positional argument
            SEXP formal_sym = TAG(f_pos);
            if (TYPEOF(formal_sym) == SYMSXP) {
                Rf_defineVar(formal_sym, arg_val, subst_env);
            }
            f_pos = CDR(f_pos);
        }
        a = CDR(a);
    }

    UNPROTECT(1); // subst_env
    return subst_env;
}

// ---------------------------------------------------------
// Simplify Pass
// ---------------------------------------------------------
// Forward declarations
SEXP Simplify_Unary(const std::string& op, SEXP arg);
SEXP Simplify_Binary(const std::string& op, SEXP left, SEXP right);

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
    const char* op = CHAR(PRINTNAME(fun));

    // Simplify all children/arguments first
    SEXP new_call = PROTECT(Rf_lcons(fun, R_NilValue));
    SEXP tail_out = new_call;
    int len = 1;

    for (SEXP tail_in = CDR(expr); tail_in != R_NilValue; tail_in = CDR(tail_in)) {
        SEXP arg_simp = PROTECT(Simplify_cpp(CAR(tail_in)));
        SETCDR(tail_out, Rf_cons(arg_simp, R_NilValue));
        tail_out = CDR(tail_out);
        
        if (TAG(tail_in) != R_NilValue) {
            SET_TAG(tail_out, TAG(tail_in));
        }
        
        UNPROTECT(1); // arg_simp
        len++;
    }

    SEXP res = new_call;

    // Control Flow Handling (`if`) on pre-simplified args
    if (std::strcmp(op, "if") == 0 && len >= 3) {
        SEXP cond = CADR(new_call);
        SEXP true_branch = CADDR(new_call);
        SEXP false_branch = (len >= 4) ? CADDDR(new_call) : R_NilValue;

        int is_true = -1;
        if (TYPEOF(cond) == LGLSXP && Rf_length(cond) > 0) {
            // logical values
            is_true = LOGICAL(cond)[0];
        } else if (TYPEOF(cond) == REALSXP && Rf_length(cond) > 0) {
            // numeric values
            is_true = (REAL(cond)[0] != 0.0);
        } else if (TYPEOF(cond) == SYMSXP) {
            // variable symbols 'T' and 'F'
            const char* sym_name = CHAR(PRINTNAME(cond));
            if (std::strcmp(sym_name, "T") == 0) is_true = 1;
            else if (std::strcmp(sym_name, "F") == 0) is_true = 0;
        }

        if (is_true == 1) {
            res = true_branch;
        } else if (is_true == 0) {
            res = false_branch; // returns R_NilValue if no else branch
        }
    } else if (len == 2) {
        res = Simplify_Unary(op, CADR(new_call));
    } else if (len == 3) {
        res = Simplify_Binary(op, CADR(new_call), CADDR(new_call));
    }

    UNPROTECT(1); // new_call
    return res;
}

// ==========================================================
// 2. UNARY RULES (Assumes 'arg' is already simplified)
// ==========================================================
SEXP Simplify_Unary(const std::string& op, SEXP arg) {
    double val;
    bool hasNum = get_numeric_val(arg, val);
    if (hasNum) {
        // whatever op, calculate the result
        SEXP call = PROTECT(Rcpp::Language(op.c_str(), arg));
        SEXP res = PROTECT(Rf_eval(call, R_BaseEnv));
        
        UNPROTECT(2); // Unprotect both call and res
        return res;
    }
    if (op == "+" || op == "(") return arg;
    
    if (op == "-") {
        // -(-y) => y (since arg is already simplified, CADR(arg) is final!)
        if (is_unary_minus(arg)) return CADR(arg); 
    }
    
    if (op == "log") {
        if (is_call_to(arg, "exp")) return CADR(arg); // y from exp(y)
    }
    if (op == "exp") {
        if (is_call_to(arg, "log") && Rf_length(arg) == 2) return CADR(arg);
    }
    
    // No rules matched: Reconstruct call
    return Rcpp::Language(op.c_str(), arg);
}

// ==========================================================
// 3. BINARY RULES (Assumes 'left' & 'right' are already simplified)
// ==========================================================
SEXP Simplify_Binary(const std::string& op, SEXP left, SEXP right) {
    double vL, vR;
    bool hasL = get_numeric_val(left, vL);
    bool hasR = get_numeric_val(right, vR);

    // Constant Folding
    if (hasL && hasR) {
        SEXP call = PROTECT(Rcpp::Language(op.c_str(), left, right));
        SEXP res = PROTECT(Rf_eval(call, R_BaseEnv));
        
        UNPROTECT(2); // Unprotect both call and res
        return res;
    }

    // Identity Rules
    if (op == "+") {
        if (hasL && vL == 0.0) return right;
        if (hasR && vR == 0.0) return left;
        
        // x + (-y) => x - y
        if (is_unary_minus(right)) {
            // SMART CONSTRUCTION: Don't call Simplify_cpp! Just call Simplify_Binary.
            return Simplify_Binary("-", left, CADR(right)); 
        }
        // -x + y => y - x
        if (is_unary_minus(left)) {
            return Simplify_Binary("-", right, CADR(left));
        }
    }
    if (op == "-") {
        if (hasL && vL == 0.0) return Simplify_Unary("-", right);
        if (hasR && vR == 0.0) return left;
        
        // x - (-y) => x + y
        if (is_unary_minus(right)) {
            return Simplify_Binary("+", left, CADR(right)); 
        }
        if (nodes_equal(left, right)) {
            return Rf_ScalarReal(0.0);
        }
    }
    
    if (op == "*") {
        if ((hasL && vL == 0.0) || (hasR && vR == 0.0)) return Rf_ScalarReal(0.0);
        if (hasL && vL == 1.0) return right;
        if (hasR && vR == 1.0) return left;
        if (hasL && vL == -1.0) return Simplify_Unary("-", right);
        if (hasR && vR == -1.0) return Simplify_Unary("-", left);
        if (hasR) return Simplify_Binary("*", right, left);
        
        // (-a) * b => -(a * b)
        if (is_unary_minus(left)) {
            SEXP a = CADR(left);
            // 1. Build a * b using the smart constructor (it folds constants if needed!)
            SEXP prod = PROTECT(Simplify_Binary("*", a, right));
            // 2. Wrap it in unary minus using the smart constructor
            SEXP res = Simplify_Unary("-", prod);
            UNPROTECT(1);
            return res;
        }
        // a * (-b) => -(a * b)
        if (is_unary_minus(right)) {
            SEXP b = CADR(right);
            // 1. Build a * b using the smart constructor (it folds constants if needed!)
            SEXP prod = PROTECT(Simplify_Binary("*", left, b));
            // 2. Wrap it in unary minus using the smart constructor
            SEXP res = Simplify_Unary("-", prod);
            UNPROTECT(1);
            return res;
        }
    }
    if (op == "/") {
        if (hasL && vL == 0.0) return Rf_ScalarReal(0.0);
        if (hasR && vR == 1.0) return left;
        if (hasR && vR == -1.0) return Simplify_Unary("-", left);
        if (hasL && vL < 0.0) return Simplify_Unary("-", Simplify_Binary("/", Rf_ScalarReal(-vL), right));
        // x / x => 1
        if (nodes_equal(left, right)) {
            return Rf_ScalarReal(1.0);
        }
        if (is_unary_minus(right)) {
            return Simplify_Unary("-", Simplify_Binary("/", left, CADR(right)));
        }
        if (is_unary_minus(left)) {
            return Simplify_Unary("-", Simplify_Binary("/", CADR(left), right));
        }
    }
    if (op == "^") {
        if (hasR) {
            if (vR == 0.0) return Rf_ScalarReal(1.0);
            if (vR == 1.0) return left;
            if (vR == 2.0) return Simplify_Binary("*", left, left);
            if (vR == 3.0) return Simplify_Binary("*", left, Simplify_Binary("*", left, left));
            if (vR < 0.) return Simplify_Binary("/", Rf_ScalarReal(1.0), Simplify_Binary("^", left, Rcpp::wrap(-vR)));
        }
    }

    if (op == "log") {
        if (nodes_equal(left, right))
            return Rf_ScalarReal(1.0);
    }

    // No rules matched: Reconstruct call
    return Rcpp::Language(op.c_str(), left, right);
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
        SEXP rules = drule.get(CHAR(PRINTNAME(fun)));
        if (rules != R_NilValue && TYPEOF(rules) == VECSXP) {
            SEXP arg_names = Rf_getAttrib(rules, R_NamesSymbol);
            int num_formals = Rf_length(rules);
            
            SEXP subst_env = PROTECT(R_NewEnv(R_EmptyEnv, 0, 0));
            
            // Look up standard formals from target environment or base package to extract default expressions
            SEXP fn_obj = env.find(op);
            //debug_print("fn_obj", fn_obj);
            // Safe formals retrieval for both closures and primitives
            SEXP fn_formals = R_NilValue;
            if (fn_obj != R_UnboundValue && fn_obj != R_NilValue) {
                if (TYPEOF(fn_obj) == CLOSXP) {
                    fn_formals = R_ClosureFormals(fn_obj);
                } else {
                    // Fallback for primitives/builtins: evaluate formals("fn_name") in base env
                    SEXP call = PROTECT(Rcpp::Language("formals", Rcpp::Language("args", op.c_str())));
                    fn_formals = PROTECT(Rf_eval(call, R_BaseEnv));
                    UNPROTECT(2);
                }
            }
            //debug_print("fn_formals", fn_formals);
            // Step 1: Pre-fill defaults from function definition
            SEXP f_ptr = fn_formals;
            while (f_ptr != R_NilValue) {
                SEXP sym = TAG(f_ptr);
                SEXP dflt = CAR(f_ptr);
                if (TYPEOF(sym) == SYMSXP && dflt != R_MissingArg) {
                    Rf_defineVar(sym, dflt, subst_env);
                }
                f_ptr = CDR(f_ptr);
            }

            // Step 2: Overlay provided actual call arguments
            SEXP curr_args = CDR(expr);
            int pos = 0;
            while (curr_args != R_NilValue) {
                SEXP arg_val = CAR(curr_args);
                SEXP arg_tag = TAG(curr_args);

                if (arg_tag != R_NilValue && TYPEOF(arg_tag) == SYMSXP) {
                    Rf_defineVar(arg_tag, arg_val, subst_env);
                } else if (pos < num_formals) {
                    SEXP formal_sym = Rf_install(CHAR(STRING_ELT(arg_names, pos)));
                    Rf_defineVar(formal_sym, arg_val, subst_env);
                }
                curr_args = CDR(curr_args);
                pos++;
            }

            // Step 3: Compute total derivative against rules
            SEXP total_deriv = R_NilValue;
            curr_args = CDR(expr);
            
            for (int i = 0; i < num_formals; i++) {
                SEXP rule_template = VECTOR_ELT(rules, i);
                // Skip arguments that have no derivative rule (e.g. lower.tail, log.p)
                if (rule_template == R_NilValue) continue;
                SEXP formal_sym = Rf_install(CHAR(STRING_ELT(arg_names, i)));
                
                // Get substituted argument or default from subst_env
                SEXP actual_arg = R_getVarEx(formal_sym, subst_env, TRUE, R_UnboundValue);;
                if (actual_arg == R_UnboundValue) continue;

                SEXP d_arg = PROTECT(deriv_cpp_internal(actual_arg, target, drule, env));
                
                if (is_numeric_val(d_arg, 0.0)) {
                    UNPROTECT(1);
                    continue;
                }
                
                SEXP df_darg = PROTECT(substitute_expr(rule_template, subst_env));
                SEXP term = PROTECT(Simplify_Binary("*", df_darg, d_arg));
                
                if (total_deriv == R_NilValue) {
                    total_deriv = term;
                    UNPROTECT(2); // df_darg, term
                } else {
                    total_deriv = PROTECT(Simplify_Binary("+", total_deriv, term));
                    UNPROTECT(3); // total_deriv_old, df_darg, term
                }
                
                UNPROTECT(1); // d_arg
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
            SEXP actual_args = CDR(expr);

            // Bind formal parameters + default values + actual call overrides
            SEXP subst_env = PROTECT(match_call_args(formals, actual_args));

            // Substitute function body using populated environment
            SEXP expanded_body = PROTECT(substitute_expr(body, subst_env));
            SEXP deriv_res = deriv_cpp_internal(expanded_body, target, drule, env);

            UNPROTECT(2); // subst_env, expanded_body
            return Simplify_cpp(deriv_res);
        }

        // --- D. UNKNOWN FUNCTION ERROR ---
        Rcpp::stop("Function '%s' has no derivative rule in 'drule' nor a function definition in 'env'", op.c_str());
    }
    //Rcout << "type top: " << Rf_type2char(TYPEOF(expr)) << std::endl;
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
