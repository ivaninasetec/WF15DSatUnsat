	!********************************************************************************************************************
	!        CLASS TO INCLUDE SOME FUNCTIONS USED IN THE INPUT CLASSESS
	!********************************************************************************************************************
	! TITLE         : 1.5D MULTILAYER FLOW
	! PROJECT       : FLOW1D HORIZONTAL SATURATED MODEL LIBRARIES
	! MODULE        : com_mod_inputs
	! URL           : ...
	! AFFILIATION   : ...
	! DATE          : ...
	! REVISION      : ... V 0.0
	! LICENSE				: This software is copyrighted 2019(C)
	!> @author
	!> Iván Campos-Guereta Díez
	!  MSc Civil Engineering by Polytechnic University of Madrid                                                     *
	!  PhD Student by University of Nottingham                                                                       *
	!  eMBA by International Institute San Telmo in Seville                                                          *
	!  ivan.camposguereta@nottingham.ac.uk
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model
	!********************************************************************************************************************

	module com_mod_inputs

	implicit none
	include 'inc_precision.fi'

	private

	public::s_com_inputs_locateblock,s_com_inputs_nextrecord
	contains

	!********************************************************************************************************************
	! S: LOCATEBLOCK(FileIndex, Bloq, Readerror)
	!--------------------------------------------------------------------------------------------------------------------
	! This subROUTine rewind the file with index "FileIndex" and search for the line with the word: 'Block '+Bloq
	! that will define the initial line for the input of that block.
	! In case of error returns Readerror<>0
	!********************************************************************************************************************

	subroutine s_com_inputs_locateblock(fileindex, bloq, readerror)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_inputs_locateblock" :: s_com_inputs_locateblock
	!DEC$ endif

	! search for 'block x' in a file and let the file pointer at that row

	integer, intent(in):: fileindex
	character*1, intent(in)::bloq
	integer, intent(out)::readerror

	character*400:: chblock
	logical::checkblock

	rewind fileindex
	readerror = 0
	checkblock=.false.
	do while (.not.checkblock)
		read (fileindex,'(A)', iostat=readerror) chblock
		!read (fileindex,'(A)') chblock
		if (readerror.ne.0) then
			write(*,*) 'BLOCK '//bloq//' not found on input file'
			stop
		end if
		if(index(chblock,'BLOCK '//bloq) .ne. 0) checkblock=.true.
	end do

	end subroutine s_com_inputs_locateblock


	subroutine s_com_inputs_nextrecord(fileindex, readerror)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_inputs_nextrecord" :: s_com_inputs_nextrecord
	!DEC$ endif

	! search for 'block x' in a file and let the file pointer at that row

	integer, intent(in):: fileindex
	integer, intent(out)::readerror

	character:: chkline
	logical::checkblock

	readerror = 0
	checkblock=.false.
	chkline='!'
	do while(chkline=='!')
		read (fileindex,'(A)', iostat=readerror) chkline
	end do
	backspace(fileindex)

	end subroutine s_com_inputs_nextrecord



	end module com_mod_inputs