from nupack import *
import numpy as np

################################################################################

def test_sample():
    a = Strand('C10G10', name='a')
    b = Strand('C10G10', name='b')
    cs = ComplexSet([a, b], complexes=SetSpec(max_size=2))
    model = Model()
    res = complex_analysis(cs, compute=['sample'], model=model, options={'num_sample': 5})
    for k, v in res.complexes.items():
        assert len(v.sample) == 5

################################################################################

def test_named():
    A1 = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A1')
    A2 = Strand('AGTCTAGGATTCGGCGTGGGTTAA', name='A2')

    assert A1 != A2 # --> False
    assert A1 == A1 # --> True

    c2 = Complex([A1, A2])
    c1 = Complex([A1, A1])

    t1 = Tube({A1: 1e-6, A2: 1e-8}, name='')
    t2 = Tube({A1: 1e-6, A2: 1e-9}, complexes=SetSpec(3, include=[c2], exclude=[c1]), name='')

################################################################################

def test_tube_analysis():
    A = Strand('CTGATCGAT', name='A')
    B = Strand('GATCGTAGTC', name='B')

    t1 = Tube({A: 1e-8, B: 1e-9}, complexes=SetSpec(3), name='0')
    t2 = Tube({A: 1e-10, B: 1e-9}, complexes=SetSpec(2), name='1')

    result = tube_analysis(tubes=[t1, t2],
        compute=['pairs', 'mfe', 'sample'], model=Model(),
        options={'num_sample': 100})

    t1 = Tube({A: 1e-8, B: 1e-8}, complexes=SetSpec(2), name='t1')
    result = complex_analysis(t1 + t2, compute=['pairs', 'mfe'], model=Model())
    result[Complex([A])] # --> ComplexResult

    t1_result = complex_concentrations(t1, result, concentrations={A: 1e-8, B: 1e-9}) # use manually specified concentrations if desired
    t2_result = complex_concentrations(t2, result) # use concentration from t2

    print(t1_result) # result concentrations

################################################################################

def test_tube_analysis_4():
    a = Strand('CAGTCGATC', name='a')
    b = Strand('ATCGACGTA', name='b')
    c = Complex([a, b])

    t1 = Tube({a: 1e-6, b: 1e-9}, complexes=SetSpec(include=[c]), name='t1')
    t2 = Tube({a: 1e-8, b: 1e-9}, complexes=SetSpec(include=[c]), name='t2')

    result = tube_analysis([t1, t2], compute=['pairs', 'mfe', 'sample', 'subopt'],
        options={'num_sample': 2, 'energy_gap': 0.5}, model=Model())
    s = str(result)
    s = result._repr_html_()

    result[t1] # --> TubeResult

    result[t1].complex_concentrations # --> [1.5e-10]
    result[t1].ensemble_pair_fractions # --> [[1.0, 0.0], [0.0, 1.0]]

    result[c]
    # pfunc for complex c
    result[c].pfunc
    # mfe for complex c
    list(result[c].mfe)
    # ppairs matrix for complex c
    result[c].pairs

################################################################################

def test_partition_function():
    a = analysis.Specification(model=Model(ensemble='some-nupack3'))
    s = RawComplex('A100')
    a.pfunc(s)
    v = a.compute()[s]

    assert v.pfunc == 1
    assert v.free_energy == 0

################################################################################

def test_ensemble_size():
    a = analysis.Specification(model=Model())
    s = RawComplex('GAAAC')
    a.ensemble_size(s)
    v = a.compute()[s]

    assert v.ensemble_size == 2

    s = RawComplex('GAAA')
    a.ensemble_size(s)
    v = a.compute()[s]

    assert v.ensemble_size == 1

################################################################################

def test_impossible_complex():
    strands = RawComplex(['A', 'A'])
    a = analysis.Specification(model=Model(ensemble='some-nupack3'))
    a.pfunc(strands)
    a.sample(strands, number=10)
    a.mfe(strands)
    a.pairs(strands)
    v = a.compute()[strands]

    assert v.pfunc == 0
    assert v.free_energy == np.inf
    assert np.all(v.pairs.to_array() == 0)
    assert len(v.mfe) == 0
    assert len(v.subopt) == 0
    assert len(v.sample) == 10
    for s in v.sample:
        assert str(s) == '.+.'
    assert v.structure_mfe() == np.inf

################################################################################

def test_matrices():
    a = analysis.Specification(model=Model(ensemble='some-nupack3'))
    s = RawComplex(['AAAAA'])
    a.pfunc(s, matrices=True)
    v = a.compute()[s]
    assert v.pfunc == 1
    assert v.pf_matrices['B'].sum() == 0

################################################################################

def test_overflow():
    s = RawComplex(['C300G300'])
    a = analysis.Specification(model=Model(ensemble='some-nupack3', material='rna95-nupack3'))
    a.pfunc(s)
    v = a.compute()[s]
    assert abs(float(v.pfunc.ln()) - 1393.4600830078125) < 0.01

    a.pairs(s)
    v = a.compute()[s]
    assert abs(v.pairs.to_array().sum(0) - 1).max() < 0.01
    assert abs(v.pairs.to_sparse().sum(0) - 1).max() < 0.01

################################################################################

