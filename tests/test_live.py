"""Pacing must slow wall time without changing the MD trajectory or saved samples."""
from types import SimpleNamespace
from unittest.mock import patch
import tempfile
import unittest
import numpy as np
from fys2160_md import Simulation, System, FCC
from fys2160_md.visualization import LiveView


class Clock:
    def __init__(self):
        self.now = 100.
        self.waits = []

    def monotonic(self):
        return self.now

    def sleep(self, delay):
        self.waits.append(delay)
        self.now += delay


class RecordingView(LiveView):
    """Replace only Jupyter transport; retain the actual pacing/snapshot methods."""
    instances = []

    def __init__(self, *, max_fps=None, frame_every=100):
        self.max_fps, self.frame_every = max_fps, frame_every
        self.widget = None
        self.last_step = None
        self.__class__.instances.append(self)

    def start(self, sim):
        from fys2160_md.visualization import time
        self.widget = SimpleNamespace(frame=self.snapshot(sim), status='running')
        self.last_step = sim.step
        self.last = time.monotonic()
        self.sent_steps = []

    def update(self, sim, **kwargs):
        previous = self.last_step
        super().update(sim, **kwargs)
        if self.last_step != previous:
            self.sent_steps.append(self.last_step)


class LiveTests(unittest.TestCase):
    def test_pacing_keeps_physics_and_recording_independent(self):
        with tempfile.TemporaryDirectory() as root:
            baseline=FCC(N=32,output_dir=root)
            a=Simulation.run(baseline,storage_name="baseline",steps=250,sample_every=70,save_every=110)
            paced=FCC(N=32,output_dir=root)
            clock=Clock()
            with patch('fys2160_md.visualization.LiveView',RecordingView), \
                 patch('fys2160_md.visualization.time.monotonic',clock.monotonic), \
                 patch('fys2160_md.visualization.time.sleep',clock.sleep):
                b=Simulation.run(paced,storage_name="paced",steps=250,sample_every=70,save_every=110,
                            show=True,max_fps=5,frame_every=100)
            view=RecordingView.instances[-1]
            self.assertEqual(view.sent_steps,[100,200,250])
            self.assertAlmostEqual(sum(clock.waits),.6)
            self.assertEqual(view.widget.status,'completed')
            self.assertEqual(view.widget.frame['step'],250)
            np.testing.assert_allclose(paced.atoms.positions,baseline.atoms.positions,atol=1e-12)
            np.testing.assert_allclose(paced.atoms.velocities,baseline.atoms.velocities,atol=1e-12)
            np.testing.assert_allclose(a.thermo,b.thermo,atol=1e-12)
            self.assertEqual([int(f['step']) for f in a.iter_frames()],
                             [int(f['step']) for f in b.iter_frames()])
            self.assertEqual(b.metadata['max_fps'],5)

    def test_full_speed_view_does_not_sleep(self):
        with tempfile.TemporaryDirectory() as root:
            sim=FCC(N=32,output_dir=root)
            clock=Clock()
            with patch('fys2160_md.visualization.LiveView',RecordingView), \
                 patch('fys2160_md.visualization.time.monotonic',clock.monotonic), \
                 patch('fys2160_md.visualization.time.sleep',clock.sleep):
                Simulation.run(sim,steps=350,show=True)
            self.assertEqual(clock.waits,[])
            self.assertEqual(RecordingView.instances[-1].widget.frame['step'],350)

    def test_interrupt_during_pacing_saves_the_completed_steps(self):
        with tempfile.TemporaryDirectory() as root:
            sim=FCC(N=32,output_dir=root)
            with patch('fys2160_md.visualization.LiveView',RecordingView), \
                 patch('fys2160_md.visualization.time.sleep',side_effect=KeyboardInterrupt):
                run=Simulation.run(sim,steps=1000,show=True,max_fps=1,frame_every=30)
            self.assertEqual(run.metadata['status'],'interrupted')
            self.assertEqual(run.metadata['completed_steps'],30)
            self.assertEqual(RecordingView.instances[-1].widget.status,'interrupted')
            self.assertEqual(RecordingView.instances[-1].widget.frame['step'],30)
            self.assertEqual(int(run.thermo.step.iloc[-1]),30)
            self.assertEqual(System.from_run(run).step,30)

    def test_invalid_display_parameters_fail_before_running(self):
        with tempfile.TemporaryDirectory() as root:
            sim=FCC(N=32,output_dir=root)
            for options in ({'max_fps':10}, {'show':True,'max_fps':0},
                            {'show':True,'max_fps':float('nan')}, {'frame_every':0}):
                with self.assertRaises(ValueError):Simulation.run(sim,**options)
            self.assertEqual(sim.step,0)


class StorageTests(unittest.TestCase):
    def test_reused_name_replaces_all_previous_output(self):
        from pathlib import Path
        with tempfile.TemporaryDirectory() as root:
            sim=FCC(N=32,output_dir=root,storage_name='my-gas')
            first=Simulation.run(sim,steps=500,save_every=50)
            self.assertEqual(len(list(first.iter_frames())),11)
            second=Simulation.run(sim,steps=10,save_every=None)
            self.assertEqual(first.path,second.path)
            self.assertEqual(second.path,Path(root).resolve()/'my-gas')
            self.assertEqual(list(second.iter_frames()),[])
            self.assertEqual(second.thermo.step.tolist(),[500,510])
            self.assertEqual(System.from_run(second).step,510)
            separate=Simulation.run(sim,steps=1,storage_name='other-gas')
            self.assertNotEqual(separate.path,second.path)
            self.assertEqual(second.thermo.step.tolist(),[500,510])

    def test_unrelated_directories_and_links_are_not_overwritten(self):
        from pathlib import Path
        with tempfile.TemporaryDirectory() as root:
            sim=FCC(N=32,output_dir=root)
            directory=Path(root)/'notes';directory.mkdir()
            (directory/'keep.txt').write_text('Important notes')
            with self.assertRaises(ValueError):Simulation.run(sim,storage_name='notes')
            self.assertEqual((directory/'keep.txt').read_text(),'Important notes')
            (Path(root)/'linked').symlink_to(directory,target_is_directory=True)
            with self.assertRaises(ValueError):Simulation.run(sim,storage_name='linked')
            for name in ['../outside','.', '/absolute', 'folder/file', 'folder\\file']:
                with self.assertRaises(ValueError):Simulation.run(sim,storage_name=name)
