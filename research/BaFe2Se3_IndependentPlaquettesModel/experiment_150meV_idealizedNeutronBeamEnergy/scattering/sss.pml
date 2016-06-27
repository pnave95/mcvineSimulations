<?xml version="1.0"?>
<!--
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! {LicenseText}
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-->

<!DOCTYPE inventory>

<inventory>

    <component name="sss">
        <facility name="sample">samples/SampleAssemblyFromXml</facility>
        <property name="dump-instrument">False</property>
        <property name="overwrite-datafiles">False</property>
        <property name="sequence">['source', 'sample', 'storage']</property>
        <property name="launcher">mpirun</property>
        <property name="output-dir">out</property>
        <facility name="storage">monitors/NeutronToStorage</facility>
        <property name="ncount">10000000.0</property>
        <property name="multiple-scattering">False</property>
        <facility name="source">sources/MonochromaticSource</facility>
        <facility name="geometer">geometer</facility>
        <property name="dump-registry">False</property>
        <property name="tracer">no-neutron-tracer</property>

        <component name="storage">
            <property name="path">scattered-neutrons</property>
        </component>


        <component name="mpirun">
            <property name="dry">False</property>
            <property name="nodelist">[]</property>
            <property name="extra"></property>
            <property name="python-mpi">`which python`</property>
            <property name="command">mpirun</property>
            <property name="debug">False</property>
            <property name="nodes-opt">-np</property>
            <property name="nodes">0</property>
        </component>


        <component name="sample">
            <property name="xml">sampleassembly/sampleassembly.xml</property>
        </component>


        <component name="source">
            <property name="probability">1.0</property>
            <property name="energy">150.0</property>
            <property name="velocity">[0.0, 0.0, 3000.0]</property>
            <property name="height">0.0</property>
            <property name="width">0.0</property>
            <property name="energy-width">0.0</property>
            <property name="time">0.0</property>
            <property name="position">[0.0, 0.0, 0.0]</property>
        </component>


        <component name="geometer">
            <property name="sample">((0, 0, 13.6), (0, 0, 0))</property>
            <property name="source">((0, 0, 13.45), (0, 0, 0))</property>
            <property name="transformer">coordinate-system-transformer</property>
            <property name="storage">((0, 0, 13.6), (0, 0, 0))</property>
            <property name="dump">False</property>
        </component>

    </component>

</inventory>

<!-- version-->
<!-- $Id$-->

<!-- Generated automatically by Renderer on Mon Jun 27 10:22:09 2016-->

<!-- End of file -->
<!-- 
 automatically created by the following command:
 $ sss -source=MonochromaticSource -dump-pml
-->

