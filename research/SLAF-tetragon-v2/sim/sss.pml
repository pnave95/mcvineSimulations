<?xml version="1.0"?>

<!DOCTYPE inventory>

<inventory>

    <component name="sss">
        <property name="sequence">['source', 'save', 'sample', 'storage']</property>
        <facility name="source">sources/NeutronFromStorage</facility>
        <facility name="save">monitors/NeutronToStorage</facility>
        <facility name="sample">samples/SampleAssemblyFromXml</facility>
        <facility name="storage">monitors/NeutronToStorage</facility>
	
        <property name="multiple-scattering">False</property>
	
        <property name="ncount">10000.0</property>

        <property name="overwrite-datafiles">False</property>
        <property name="output-dir">out</property>
	

        <component name="source">
            <property name="path">beam/out/neutrons</property>
        </component>


        <component name="save">
            <property name="path">saved-neutrons</property>
        </component>


        <component name="sample">
 	    <property name="xml">sampleassembly/sampleassembly.xml</property>
        </component>


        <component name="storage">
            <property name="path">scattered-neutrons</property>
        </component>


        <component name="geometer">
            <property name="source">((0, 0, -0.15), (0, 0, 0))</property>
            <property name="save">((0, 0, -0.149), (0, 0, 0))</property>
            <property name="sample">((0, 0, 0), (0, 0, 0))</property>
            <property name="storage">((0, 0, 0), (0, 0, 0))</property>
        </component>


    </component>

</inventory>

<!-- End of file -->

