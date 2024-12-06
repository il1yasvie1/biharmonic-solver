from firedrake import *


class BiharmonicProblem:
    def __init__(self, cfg):
        self.config = cfg
        self.method = cfg['method']
        self.mesh = self._create_mesh()
        self.function_space = self._create_function_space()

    def _create_mesh(self):
        supported_mesh_types = {
            'RectangleMesh': RectangleMesh,
            'SquareMesh': SquareMesh,
        }
        mesh_cfg = self.config['mesh']
        mesh_type = mesh_cfg['type']
        if mesh_type in supported_mesh_types:
            mesh = supported_mesh_types[mesh_type](**mesh_cfg['parameters'])
            return mesh
        else:
            raise NotImplementedError

    def _create_function_space(self):
        fs_cfg = self.config['function_space']
        if self.method == 'c0':
            fs = VectorFunctionSpace(self.mesh, 'CG', degree=fs_cfg['V']['degree'], dim=3)

        elif self.method == 'dg':
            fs = VectorFunctionSpace(self.mesh, 'DG', degree=fs_cfg['V']['degree'], dim=3)

        elif self.method == 'mix':
            V = VectorFunctionSpace(self.mesh, 'CG', degree=fs_cfg['V']['degree'], dim=3)
            Q = TensorFunctionSpace(self.mesh, 'CG', degree=fs_cfg['Q']['degree'], shape=(3, 2))
            fs = V*Q*Q

        else:
            raise NotImplementedError
        return fs

    def solve(self):
        bcs_cfg = self.config['boundary_conditions']
        mesh = self.mesh
        x = SpatialCoordinate(mesh)
        n = FacetNormal(mesh)
        h = avg(CellVolume(mesh)) / FacetArea(mesh)

        f_expr = as_vector([eval(expr) for expr in self.config['f']])
        g_expr = as_vector([eval(expr) for expr in bcs_cfg['g']])

        if bcs_cfg['phi']:
            phi_expr = as_matrix([[eval(expr) for expr in row] for row in bcs_cfg['phi']])

        sub_domain = tuple(bcs_cfg['sub_domain'])
        if self.config['solver_parameters']:
            solver_parameters = dict(self.config['solver_parameters'])
        else:
            solver_parameters = {}

        r = Constant(bcs_cfg['r'])

        if self.method in {'c0', 'dg'}:
            V = self.function_space
            f = Function(V).interpolate(f_expr)
            g = Function(V).interpolate(g_expr)
            if bcs_cfg['phi']:
                phi = phi_expr
            y = Function(V)

            if self.config['initial_guess']:
                y0_expr = as_vector([eval(expr) for expr in self.config['initial_guess']['y0']])
                y0 = Function(V).interpolate(y0_expr)
                y.assign(y0)

            E = (inner(grad(grad(y)), grad(grad(y))) - dot(y, f))*dx
            E += dot(dot(jump(grad(y)), n('+')), dot(dot(avg(grad(grad(y))), n('+')), n('+')))*dS
            E += r*dot(dot(jump(grad(y)), n('+')), dot(jump(grad(y)), n('+')))/(2*h)*dS

            if bcs_cfg['phi']:
                E += dot(dot(grad(y)-phi, n), dot(dot(grad(grad(y)), n), n))*ds(sub_domain)
                E += r*dot(dot(grad(y)-phi, n), dot(grad(y)-phi, n))/(2*h)*ds(sub_domain)

            if self.method == 'dg':
                r0 = Constant(bcs_cfg['r0'])
                E += dot(dot(avg(div(grad(grad(y)))), n('+')), jump(y))*dS
                E += dot(dot(div(grad(grad(y))), n), y-g)*ds(sub_domain)
                E += r0*dot(jump(y), jump(y))/(2*h**(3/2))*dS
                E += r0*dot(y-g, y-g)/(2*h**(3/2))*ds(sub_domain)
                F = derivative(E, y)
                solve(F == 0, y, solver_parameters=solver_parameters)

            elif self.method == 'c0':
                F = derivative(E, y)
                bc = DirichletBC(V, g, sub_domain)
                solve(F == 0, y, bcs=[bc], solver_parameters=solver_parameters)

            return y

        elif self.method == 'mix':
            Z = self.function_space
            f = Function(Z.sub(0)).interpolate(f_expr)
            g = Function(Z.sub(0)).interpolate(g_expr)
            if bcs_cfg['phi']:
                phi = Function(Z.sub(1)).interpolate(phi_expr)
            z = Function(Z)

            if self.config['initial_guess']:
                z0 = Function(Z)
                y0, u0, p0 = z0.subfunctions
                if self.config['initial_guess']['y0']:
                    y0_expr = as_vector([eval(expr) for expr in self.config['initial_guess']['y0']])
                    y0 = Function(Z.sub(0)).interpolate(y0_expr)
                if self.config['initial_guess']['u0']:
                    u0_expr = as_matrix([[eval(expr) for expr in row] for row in self.config['initial_guess']['u0']])
                    u0 = Function(Z.sub(1)).interpolate(u0_expr)
                if self.config['initial_guess']['p0']:
                    p0_expr = as_matrix([[eval(expr) for expr in row] for row in self.config['initial_guess']['p0']])
                    p0 = Function(Z.sub(2)).interpolate(p0_expr)
                z.assign(z0)

            y, u, p = split(z)
            E = (inner(grad(u), grad(u)) + inner(grad(grad(y)), grad(grad(y))))*dx
            E -= dot(y, f)*dx
            E += inner(p, grad(y) - u)*dx
            E -= inner(dot(avg(grad(u)), n('+')), jump(u))*dS
            E += r*inner(jump(u), jump(u))/(2*h)*dS

            if bcs_cfg['phi']:
                E -= inner(dot(grad(u), n), u - phi)*ds(sub_domain)
                E += r*inner(u - phi, u - phi)/(2*h)*ds(sub_domain)
            F = derivative(E, z)

            bc = DirichletBC(Z.sub(0), g, sub_domain)
            solve(F==0, z, bcs=[bc], solver_parameters=solver_parameters)
            return z.sub(0)
        else:
            raise NotImplementedError
